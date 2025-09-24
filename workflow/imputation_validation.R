# 🟡 Library ----
library(tidyverse)
library(tidytext)
library(ggcorrplot)
library(corrr)
library(caret)
library(rcompanion)
library(magrittr)
library(terra)
library(smoothr)
library(exactextractr)
library(sf)
library(pfif)
library(pals)
library(randomForest)





# 🟡 Functions ----
rstudioapi::getSourceEditorContext()$path %>% 
  dirname() %>% 
  gsub("/eFRI_workflow", ., replacement = "") %>% 
  paste0("/functions/imputation") %>% 
  list.files(full.names = T) %>% 
  map(source)





# 🟡 Setwd and parameters ----
sector <- "RMF"
setwd(paste0("D:/00_Ontario_eFRI/", sector))
rm(sector)
dem <- rast("./metrics/lidar/dem.tif")
epsg <- st_crs(dem)





# 🟡 List of data ----
rstudioapi::getSourceEditorContext()$path %>% 
  dirname() %>% 
  gsub("/eFRI_workflow", ., replacement = "") %>% 
  paste0("/functions/other/read_data.R") %>% 
  source()





# 🟡 Read base data ----
# Catalog
st_read("./ctg/ctg.shp") %>% 
  st_transform(epsg) -> ctg

# Forest inventory polygons (fri)
# set.seed(1) # SÉLECTION ALÉATOIRE DE 10 % DES ENTITÉES

st_read("./shp/PolygonForest.shp") %>% 
  rowid_to_column("id") %>% 
#  sample_frac(0.005) %>% # SÉLECTION ALÉATOIRE DE 1 % DES ENTITÉES
  st_transform(epsg) -> fri_polygons

# All metrics in a spat raster with mutiple layers
metrics_infos %>%
  filter(type != "dendrometric",
         res == 20) %T>% 
  {pull(.,name) ->> metrics_names} %>% # Extract metrics names in the right order
  pull(path) %>% 
  map(rast) %>% 
  map(project, dem, method = "bilinear") %>% 
  map(resample, dem) %>% 
  rast -> metrics

# Assign metrics names
names(metrics) <- metrics_names

# Landcover
rast("./metrics/other/landcover.tif") -> landcover

# Forest_Fire_1985-2020
rast("./metrics/other/forest_fire_1985_2020.tif") %>% 
  as.polygons() %>% 
  st_as_sf() %>% 
  dplyr::filter(.[[1]] != 0) %>% 
  st_cast("MULTIPOLYGON") %>% 
  st_cast("POLYGON") %>% 
  rename(year_fire = 1) -> forest_fire_1985_2020_poly


# Forest_Harvest_1985-2020
rast("./metrics/other/forest_harvest_1985_2020.tif") %>% 
  as.polygons() %>% 
  st_as_sf() %>% 
  dplyr::filter(.[[1]] != 0) %>% 
  st_cast("MULTIPOLYGON") %>% 
  st_cast("POLYGON") %>% 
  rename(year_harvest = 1) -> forest_harvest_1985_2020_poly

rm(dem, 
   epsg, 
   metrics_infos)





# 🟡 FRI polygons ----
# Find ids that had an harvest between inventory and lidar flight so that the inventory is not good
fri_polygons %>%
  dplyr::select(id, 
                year = YRUPD) %>% 
  st_intersection(ctg %>% dplyr::select(year_lidar = Fl_Cr_Y)) %>% 
  st_intersection(forest_harvest_1985_2020_poly %>% dplyr::select(year_harvest)) %>% 
  mutate(drop = ifelse(year_harvest > year & year_harvest < year_lidar, 1, 0)) %>% 
  filter(drop == 1) %>% 
  pull(id) %>% 
  unique -> drop_id_harvest

# Find ids that had a fire between inventory and lidar flight so that the inventory is not good
fri_polygons %>%
  dplyr::select(id, 
                year = YRUPD) %>% 
  st_intersection(ctg %>% dplyr::select(year_lidar = Fl_Cr_Y)) %>% 
  st_intersection(forest_fire_1985_2020_poly %>% dplyr::select(year_fire)) %>% 
  mutate(drop = ifelse(year_fire > year & year_fire < year_lidar, 1, 0)) %>% 
  filter(drop == 1) %>% 
  pull(id) %>% 
  unique -> drop_id_fire

fri_polygons %<>% 
  filter(!id %in% c(drop_id_harvest, drop_id_fire)) %>% 
  mutate.metrics(metrics, fun = "median") %>%   
  mutate.landcover(landcover) %>%   
  mutate.prop.forested() %>% 
  filter(POLYTYPE == "FOR") %>% # remove polygon that are not forested
  filter(PROPFORESTED >= 50) %>% # remove polygon that have at least 50 % of forested area
  filter(z_p95 >= 5, z_above2 >= 50) %>%  # remove polygon with only low vegetation according to lidar data
  mutate.species.prop("SPCOMP") %>% # RMF
  mutate.species.order("SPCOMP") %>% # RMF
  mutate.forest.type %>% 
  mutate.functionals.groups

# Add coordinates of centroid
fri_polygons %<>% 
  st_geometry() %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  bind_cols(fri_polygons, .)

rm(ctg,
   landcover,
   forest_fire_1985_2020_poly,
   forest_harvest_1985_2020_poly,
   drop_id_harvest,
   drop_id_fire,
   metrics)





# 🟡 Correlation matrix and variable selection ----
fri_polygons %>% 
  st_drop_geometry() %>% 
  dplyr::select(-id) %>% 
  dplyr::select(any_of(metrics_names)) %>%
  cor(use = "complete.obs") -> corr_matrix

metrics_names[!metrics_names %in% findCorrelation(corr_matrix, cutoff = 0.80, names = TRUE)] -> metrics_names_filtered

# Remove irrelevents metrics
metrics_names_filtered[!metrics_names_filtered %in% c("EVI", "NDMI", "twi")] -> metrics_names_filtered

# Add relevent missing metrics
c(metrics_names_filtered,
  c("z_p50",
    "z_p80",
    "z_mean",
    "z_sd",
    "B2",
    "B3",
    "B5",
    "NDVI")) -> metrics_names_filtered

as.dist(1 - corr_matrix) %>% # Calculate distance matrix (1 - correlation)
  hclust(method = "ward.D") %>% # Perform hierarchical clustering
  {.$order} %>% # Get the order of the clusters
  {corr_matrix[., .]} %T>% # Reorder the matrix
  {rownames(.) ->> corr_matrix_levels} %>% # Keep in global environment the levels to maintain clustering order in the plot
  as.data.frame() %>%
  rownames_to_column("Variable1") %>%
  pivot_longer(-Variable1, # Convert the reordered matrix to long format for plotting
               names_to = "Variable2", 
               values_to = "Correlation") %>% 
  mutate(Variable1 = factor(Variable1, levels = corr_matrix_levels), # Set factor levels to maintain clustering order in the plot
         Variable2 = factor(Variable2, levels = corr_matrix_levels)) -> corr_matrix_long

# Create the heatmap plot with ggplot2
ggplot(data = corr_matrix_long, 
            mapping = aes(x = Variable1, y = Variable2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_tile(data = corr_matrix_long %>% filter(Variable1 %in% metrics_names_filtered & Variable2 %in% metrics_names_filtered), 
            mapping =  aes(Variable2, Variable1), color = "black", size = 2) + # Select filtered variables and display them with black rectangles
  geom_text(aes(label = round(Correlation, 2)), 
            color = "black", size = 3) +
  geom_rect(size=1, fill=NA, colour="black",
            aes(xmin=0.5, xmax=0.5, ymin=0.5, ymax=0.5)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Reordered Correlation Matrix (Hierarchical Clustering)", x = "", y = "")

ggsave("./results/correlation_matrix.jpg",
       width=20, 
       height=20, 
       dpi=300)

rm(corr_matrix,
   corr_matrix_levels,
   corr_matrix_long)





# 🟡 PCA ----
rasters_tibble %>%
  sample_frac(0.01) %>% 
  drop_na %>%
  prcomp(scale = TRUE) -> pca
summary(pca)
pca$rotation

ggbiplot(pca, 
         obs.scale = 1, 
         var.factor = 4,
         var.axes = TRUE,
         varname.color = metrics_infos %>% filter(name %in% (rasters_tibble %>% names)) %>% pull(color)) 

ggsave("./figures/pca_trim.pdf", 
       width=20, 
       height=20, 
       dpi=300)





# 🟡 Imputation on target variables, on knn variable at the time ----
# Building the combinaisons
c(combn(metrics_names_filtered, 1, simplify = F),
  combn(metrics_names_filtered, 2, simplify = F),
  combn(metrics_names_filtered, 3, simplify = F),
  combn(metrics_names_filtered, 4, simplify = F),
  combn(metrics_names_filtered, 5, simplify = F)) -> knn_vars_combs

# Find correlated groups
fri_polygons %>% 
  st_drop_geometry() %>% 
  dplyr::select(-id) %>% 
  dplyr::select(any_of(metrics_names_filtered)) %>%
  cor(use = "complete.obs") -> corr_matrix

corr_groups <- creditmodel::get_correlation_group(corr_matrix, p = 0.6) %>%
  keep(~ length(.x) > 1) # Keep groups that have more than 1 variable

# Filter combinations - keep only those without correlation
# ".x" refer to one combinaison and "group" to one correlated group. For each combinaison (knn_vars_combs),
# each correlated group (corr_groups) are analysed if the combinaison contain less or  
# equal 1 metrics of the group (sum() <= 1). A logical vector is returned and if 
# all are TRUE (all()), the combinaison is keep (keep()) if the functions didn't
# find in the combinaison 2 metrics of more of the same group.
knn_vars_combs_filtered <- knn_vars_combs %>%
  keep(~ all(map_int(corr_groups, function(group) sum(.x %in% group)) <= 1))

rm(corr_matrix,
   corr_groups)

results <- list()

r <- 1

for(target_variable_selected in c("SP_NO_1", "SP_NO_2", "FUNCTIONAL_GROUP_3", "FUNCTIONAL_GROUP_5")){
  for(knn_vars in knn_vars_combs_filtered){
    # knn_vars <- c("z_mean", "rumple_index", "z_p80", "z_sd", "B6") # Ethan comb (pcum_80 in original and not z_p80)
    # knn_vars <- c("rumple_index","sagawi","z_p95","B6") # Custom comb
    
    # Add X and Y
    
    knn_vars <- c(knn_vars,
                  "X",
                  "Y")
    
    cat(paste0("Target variable : ", target_variable_selected, " / Knn variables : ", paste(knn_vars, collapse = " - "), "\n"))
    
    # Select only target columns of the fri polygons with values (drop_na)
    fri_polygons %>%
      dplyr::select(any_of(c("id", 
                             knn_vars, 
                             target_variable_selected))) %>% # The unique id, all variables used for knn imputation and the current target variable selected
      mutate(!!sym(target_variable_selected) := factor(!!sym(target_variable_selected))) %>% 
      drop_na(everything()) -> fri_polygons_temp
    
    # Perform knn imputation
    knn.inputation(reference_polygons = fri_polygons_temp, # With all fri polygons with data (non-na)
                   target_polygons = fri_polygons_temp, # For all grm polygons with data (non-na)
                   knn_variables = knn_vars, # With the 5 most correlated variable
                   target_variables = target_variable_selected, # For the current target variable selected
                   k = 5) -> imputation_results
    
    # Get scores
    confusionMatrix(data = fri_polygons_temp %>% pull(target_variable_selected), 
                    reference = imputation_results %>% pull(target_variable_selected)) -> scores
    
    # Macro-average:
    #   
    #   Averages metric equally across all classes.
    # 
    # Useful when you care equally about all classes regardless of support.
    # 
    # Micro-average:
    #   
    #   Aggregates all TP, FP, FN across classes and then computes metrics.
    # 
    # Useful when you care more about overall performance, especially with imbalanced classes.
    #
    # For multi-class classification, micro-average precision equals micro-average recall and equals accuracy.
    
    scores[["byClass"]] %>% 
      as_tibble %>% 
      replace(is.na(.), 0) %>% 
      mutate(true_positive = scores[["table"]] %>% diag(),
             false_positive = rowSums(scores[["table"]]) - true_positive,
             false_negative = colSums(scores[["table"]]) - true_positive) -> scores_by_class
    
    tibble(target_var = target_variable_selected,
           knn_vars = paste(knn_vars, collapse = ","),
           accuracy = scores[["overall"]][["Accuracy"]],
           kappa = scores[["overall"]][["Kappa"]],
           macro_precision = scores_by_class$Precision %>% mean,
           macro_recall = scores_by_class$Recall %>% mean,
           macro_F1 = scores_by_class$F1 %>% mean) -> results_temp # A VALIDER AVEC ALEXANDRE SI ON PEUT METTRE 0 DANS RECALL ET FAIRE UNE MOYENNE DES F1
    
    results_temp
    
    results[[r]] <- results_temp
    
    r <- r + 1
  }
}

results %<>% 
  bind_rows()

results %>% 
  saveRDS("./results/imputation_results.rds")

rm(r,
   target_variable_selected,
   knn_vars,
   fri_polygons_temp,
   imputation_results,
   scores,
   scores_by_class,
   results_temp)


# 🟡 Display results ----
# Results most frequent variables in top 20 models
results %>% 
  group_by(target_var) %>% 
  slice_max(accuracy, n = 20) %>%
  separate_rows(knn_vars, sep = ",") %>%
  filter(!knn_vars %in% c("X", "Y")) %>% 
  count(knn_var = knn_vars, sort = TRUE) %>% 
  arrange(target_var) -> results_n_var

results_n_var %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(knn_var, n, target_var), y = n), fill = "steelblue") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~target_var, scales = "free_y") +
  labs(title = "Most frequent variables in top 20 models",
       x = "Variable", y = "Occurence") +
  theme_minimal()

ggsave("./results/most_frequent_variables.jpg",
       width=20, 
       height=10, 
       dpi=300)

# Results mean performance
results %>% 
  separate_rows(knn_vars, sep = ",") %>%  
  filter(!knn_vars %in% c("X", "Y")) %>% 
  rename(knn_var = knn_vars) %>% 
  group_by(knn_var,
           target_var) %>% 
  summarise(accuracy = mean(accuracy),
            kappa = mean(kappa),
            macro_precision = mean(macro_precision),
            macro_recall = mean(macro_recall),
            macro_F1 = mean(macro_F1)) %>% 
  arrange(target_var,
          desc(accuracy)) -> results_performance

results_performance %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(knn_var, accuracy, target_var), y = accuracy), fill = "steelblue") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~target_var, scales = "free_y") +
  labs(title = "Mean accuracy by variable",
       x = "Variable", y = "Mean accuracy") +
  theme_minimal()

ggsave("./results/mean_accuracy.jpg",
       width=20, 
       height=10, 
       dpi=300)

results_performance %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(knn_var, kappa, target_var), y = kappa), fill = "steelblue") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~target_var, scales = "free_y") +
  labs(title = "Mean kappa by variable",
       x = "Variable", y = "Mean kappa") +
  theme_minimal()

ggsave("./results/mean_kappa.jpg",
       width=20, 
       height=10, 
       dpi=300)

results_performance %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(knn_var, macro_precision, target_var), y = macro_precision), fill = "steelblue") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~target_var, scales = "free_y") +
  labs(title = "Mean macro precision by variable",
       x = "Variable", y = "Mean macro precision") +
  theme_minimal()

ggsave("./results/mean_macro_precision.jpg",
       width=20, 
       height=10, 
       dpi=300)

results_performance %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(knn_var, macro_recall, target_var), y = macro_recall), fill = "steelblue") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~target_var, scales = "free_y") +
  labs(title = "Mean macro recall by variable",
       x = "Variable", y = "Mean macro recall") +
  theme_minimal()

ggsave("./results/mean_macro_recall.jpg",
       width=20, 
       height=10, 
       dpi=300)

results_performance %>% 
  ggplot() +
  geom_col(aes(x = reorder_within(knn_var, macro_F1, target_var), y = macro_F1), fill = "steelblue") +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~target_var, scales = "free_y") +
  labs(title = "Mean macro F1 by variable",
       x = "Variable", y = "Mean macro F1") +
  theme_minimal()

ggsave("./results/mean_macro_F1.jpg",
       width=20, 
       height=10, 
       dpi=300)

# Results best models
results %>% 
  group_by(target_var) %>% 
  slice_max(accuracy, n = 1) -> results_best_model

results_best_model %>% 
  write.csv("./results/best_model.csv")

# Results with a random forest
# Extract the 5 variables that have the more correlation with the current target variable selected
results %>%
  rowid_to_column("comb_id") %>% 
  mutate(knn_var = strsplit(knn_vars, ",")) %>%
  unnest(knn_var) %>% 
  mutate(value = 1) %>%
  pivot_wider(names_from = knn_var, # Pivot wider to get 0/1 for presence/absence
              values_from = value,
              values_fill = 0) %>% 
  filter(target_var == "FUNCTIONAL_GROUP_3") %>% 
  dplyr::select(any_of(c(metrics_names_filtered, "accuracy"))) -> results_tri

rf_model <- randomForest(accuracy ~ ., data = results_tri, importance = TRUE)

rf_model %>% 
  importance(type = 1) %>% # Variable importance
  as.data.frame() %>% 
  rownames_to_column("variable") %>% 
  ggplot(aes(x = reorder(variable, `%IncMSE`), y = `%IncMSE`)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Variable Importance (Random Forest)",
       x = "Variable", y = "Increase in MSE (%)") +
  theme_minimal()

ggsave("./results/random_forest_variable_importance.jpg",
       width=20, 
       height=10, 
       dpi=300)

rm(results_n_var,
   results_performance,
   results_best_model,
   results_tri,
   rf_model)





# # SUREMENT A SUPPRIMER
# # RF directly on polygon
# fri_polygons %>% 
#   st_drop_geometry() %>% 
#   dplyr::select(any_of(c(metrics_names_filtered, 
#                          target_variable_selected))) %>% 
#   mutate(!!sym(target_variable_selected) := factor(!!sym(target_variable_selected))) %>% 
#   drop_na(everything()) -> fri_polygons_selected
# 
# rf <- randomForest(FUNCTIONAL_GROUP_3 ~ ., data = fri_polygons_selected, importance = TRUE)