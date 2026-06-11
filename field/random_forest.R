# Questions
# Pourquoi on fait une corrélation cutoff de 0.9 pour connaitre les variables a éléminer, mais au final on garde tout ?
# Pourquoi on fait deux RF pour le nombre optimal d'arbres ?
# Pourquoi on fait le nombre optimal d'arbres alors qu'il existe une fonction intégrée qui le fait plus rapidement ? (opt_importance())
# Pourquoi 7000 arbres alors que la fonction et l'analyse de sensibilité en propose 5000 ?
# Pourquoi ne pas mettre le set.seed juste avant le train ou bien dans le ctrl du train ?
# Pourquoi les arbres ne sont pas faits dans ranger ?
# Pourquoi faire un k-cross-fold ? -> tout le monde fait ça
# Est-ce que la permutation a été envisagée pour retirer itérativement les variables du modèles ? La simple importance ne fonctionne pas lorsqu'il existe plusieurs variable corrélés entre elles
# Pourquoi le fractional cover n'a pas été utilisé ?
# Faire le r de pearson avec la variable a prédire pour choisir la variable dans les clusters

# PCA avec corrélation avec la variable réponse (prcomp)

# Johanne White 2017 - Propose une méthode?
# For modelling of most forest inventory
# attributes with good predictive performance, a few metrics
# representing plot height, variability of plot height, and plot
# canopy cover are often reliable predictors for basal area,
# volume and biomass (Table 8; Lefsky et al. 2005; Li et al.
#                     2008; Bouvier et al. 2015).
#
# The aim of model development, regardless
# of method, should be parsimony (e.g. selecting only a
#                                small set of relevant predictors).
#
# To reiterate, in general, metrics informing
# on height, the variability in height, and the amount of
# vegetation present (as a minimum) should support a range
# of applications and can serve as a useful starting point for
# model development.




# Librarie, functions and parameters ----
library(tidyverse)
library(randomForest)
library(plotly)
library(CAST)
library(caret)
library(doParallel)
library(optRF)
library(pheatmap)
library(terra)
library(vip)

filter_lidar_parameters <- function(data,
                                    radius = "1128",
                                    zmin = "NA"){
  data %>%
    select(Name, contains(paste0("_radius_", radius, "_cm"))) %>%
    rename_with(.fn = ~ gsub(paste0("_radius_", radius, "_cm"), "", .)) %>%
    select(Name, contains(paste0("_zmin_", zmin, "_cm"))) %>%
    rename_with(.fn = ~ gsub(paste0("_zmin_", zmin, "_cm"), "", .)) -> data_filtered

  return(data_filtered)
}

filter_ufc_parameters <- function(data,
                                  radius = "1128"){
  data %>%
    select(Name, contains(paste0("_radius_", radius))) %>%
    rename_with(.fn = ~ gsub(paste0("_radius_", radius), "", .)) -> data_filtered

  return(data_filtered)
}

correct_placettes_names <- function(data){

  data %>%
    mutate(Name = gsub("_centre", "", Name)) %>% # enleve l'appelation centre
    mutate(Name = gsub("X", "", Name)) %>% # enleve les X dans les noms de parcelles qui indiquent que la parcelle a due etre deplacee
    filter(!grepl('FC', Name)) %>% # enleve la parcelle FC
    mutate(Name = case_when(Name == 'S10_52' ~ 'S10_5', .default = Name)) -> data_corrected # corrige un nom de parcelle

  return(data_corrected)
}

rename_pz <- function(data){

  data %>%
    rename("pz_0-0.15" = "pz_0.0.15",
           "pz_0.15-2" = "pz_0.15.2",
           "pz_2-5" = "pz_2.5",
           "pz_5-10" = "pz_5.10",
           "pz_10-20" = "pz_10.20",
           "pz_20-30" = "pz_20.30") -> data_renamed

  return(data_renamed)

}

# Les plus complexes qui ont ete mise de cotes pour faciliter le calcul des metriques
selected_var <- c("zcv",
                  "zskew",
                  "zq1",
                  "zq5",
                  "zq10",
                  "zq15",
                  "zq20",
                  "zq25",
                  "zq30",
                  "zq35",
                  "zq40",
                  "zq45",
                  "zq50",
                  "zq55",
                  "zq60",
                  "zq65",
                  "zq70",
                  "zq75",
                  "zq80",
                  "zq85",
                  "zq90",
                  "zq95",
                  "zq99",
                  "pzabovemean",
                  "pzabove2",
                  "pzabove5",
                  "zpcum1",
                  "zpcum2",
                  "zpcum3",
                  "zpcum4",
                  "zpcum5",
                  "zpcum6",
                  "zpcum7",
                  "zpcum8",
                  "zpcum9",
                  "pz_below_0",
                  "pz_0-0.15",
                  "pz_0.15-2",
                  "pz_2-5",
                  "pz_5-10",
                  "pz_10-20",
                  "pz_20-30",
                  # "pz_above_30", # pz_above_30 elimine car ne contient que des 0 pour RMF
                  "ziqr",
                  "zMADmean",
                  "zMADmedian",
                  "CRR",
                  # "zentropy", # force d'enlever zentropy car valeur manquante pour une des parcelles improductives de OVF
                  "ufc_circle_1m_to_2m",
                  "ufc_circle_2m_to_3m",
                  "ufc_circle_3m_to_4m",
                  "ufc_circle_1m_to_4m",
                  "aspect",
                  "dem",
                  "relative_topographic_position_60m",
                  "relative_topographic_position_100m",
                  "relative_topographic_position_200m",
                  "relative_topographic_position_400m",
                  "relative_topographic_position_800m",
                  "sagawi",
                  "slope",
                  "twi")

# Random forest parameters
# Hyperparameter tuning grid
tuneGrid_rf <- expand.grid(.mtry = c(2, 5, 10, 15)) #Reçoit avertissement quand le nombre de variables est en-dessous des valeurs de mtry mais après vérification ça ne semble pas significativement affecter les résultats

pred_wrapper <- function(object, newdata) {predict(object, newdata, type = "raw")}

threshold_corr <- 0.9

sites <- c("RMF", "OVF")

variables <- c("vmb_ha", "st", "dhpq", "dens")





# Setwd ----
setwd("D:/00_Ontario_eFRI/random_forest")





# Importation et préparation des données ----
# Donnees topo
readRDS("./base_data/placettes_data_improductives_topo.rds") %>%
  select(Name,
         aspect,
         dem,
         relative_topographic_position_100m,
         relative_topographic_position_200m,
         relative_topographic_position_400m,
         relative_topographic_position_800m,
         relative_topographic_position_60m,
         sagawi,
         slope,
         twi) -> imp_topo

readRDS("./base_data/placettes_data_topo.rds") %>%
  correct_placettes_names() %>%
  bind_rows(imp_topo) -> topo

# Metriques lidar ufc
readRDS("./base_data/placettes_data_improductives_topo.rds") %>%
  filter_ufc_parameters() %>%
  dplyr::select(Name, contains("ufc")) -> imp_ufc

readRDS("./base_data/placettes_data.rds") %>%
  filter_ufc_parameters() %>%
  dplyr::select(Name, contains("ufc")) %>%
  correct_placettes_names() %>%
  bind_rows(imp_ufc) -> ufc


# Metriques lidar set3
readRDS("./base_data/placettes_data_improductives_topo.rds") %>%
  filter_lidar_parameters() %>%
  separate_wider_delim(Name, delim = "_", names = c("Bloc", "Plot"), cols_remove = FALSE) -> imp_lidar

readRDS("./base_data/placettes_data.rds") %>%
  filter_lidar_parameters() %>%
  rename_pz() %>%
  correct_placettes_names() %>%
  separate_wider_delim(Name, delim = "_", names = c("Bloc", "Plot"), cols_remove = FALSE) %>%
  bind_rows(imp_lidar) -> lidar

lidar %>%
  full_join(ufc, by = "Name") %>%
  full_join(topo, by = "Name") -> data

rm(imp_topo,
   topo,
   imp_ufc,
   ufc,
   imp_lidar,
   lidar)

# Dendrometriques
read.csv("./base_data/dendro_ontario.csv", sep = ';') %>%
  dplyr::select(Name = pe,
                dens = Densite_tiges_ha_Marchandes_Vivante,
                st = Surface_terr_ha_Marchandes_Vivante,
                dhpq = DHQ_cm_Marchandes_Vivante,
                vmb) %>%
  filter(Name != "S7_4") %>% #enlève une ligne avec un nom de parcelle probablement erroné (sans match dans les metriques)
  mutate(vmb_ha = vmb * 25) %>%
  full_join(data, by = "Name") %>%
  mutate(Site = case_when(Bloc %in% paste0("S", 1:5) ~ "RMF",
                          Bloc %in% paste0("S", 6:10) ~ "OVF")) -> data





# # Visualisation des donnees ----
# data %>%
#   ggplot() +
#   geom_point(aes(rumple, vmb_ha, color = Bloc)) +
#   facet_wrap(~Site)
#
#
#

# # Correlation entre les variables pour la selection PY ----
# # s <- "RMF"
# for(s in sites){
#   print(s)
#   data %>%
#     filter(Site == s) %>%
#     dplyr::select(all_of(selected_var)) %>%
#     cor() %>%
#     findCorrelation(cutoff = 0.9, name = TRUE) %>%
#     print()
# }
#
#
#
#

# # Correlation entre les variables via une pca pour la selection AMB ----
# # s <- "OVF"
# for(s in sites){
#   print(s)
#
#   data %>%
#     filter(Site == s) %>%
#     dplyr::select(all_of(selected_var)) -> data_pca
#
#   data_pca[complete.cases(data_pca), , drop = FALSE] -> data_pca
#
#   data_pca %>%
#     stats::prcomp(center = TRUE,
#                   scale. = TRUE) -> pca
#
#   pca$sdev^2 / sum(pca$sdev^2) -> var_explained
#
#   data.frame(PC = paste0("PC", seq_along(var_explained)),
#              Variance_Explained = var_explained * 100,
#              Cumulative = cumsum(var_explained) * 100) -> var_table
#
#   k <- sum(var_table$Cumulative < 90)   # je suggérerais de retenir le nombre de composantes qui permet d'expliquer 90 ou 80% de la variance totale
#
#   # Methode FL
#   var_exp <- pca$sdev^2 / sum(pca$sdev^2)
#   importance <- rowSums(pca$rotation[, 1:k, drop = FALSE]^2 * var_exp[1:k])
#   sort(importance, decreasing = TRUE)
#
#   # Rotation varimax
#   rotation <- stats::varimax(pca$rotation[, 1:k, drop = FALSE])
#   rot_loadings <- as.matrix(rotation$loadings)
#
#   # Visualisation des loadings
#   pheatmap::pheatmap(rot_loadings,
#                      main = "Rotated PCA loadings")
#
#   # Sélection des variables les plus représentatives par composante
#   n_var <- 3
#   # x <- 3
#
#   map(seq_len(k), # Nombre de composante à regarder
#           function(x){
#
#             rot_loadings[, x] %>%
#               abs() %>%
#               sort(decreasing = TRUE) %>%
#               names() %>%
#               {.[seq_len(n_var)]}
#
#           }) %>%
#     unlist() %>%
#     unique() %>%
#     print()
# }
#
#
#
#

# # Test nombre optimal d'arbres pour ovf et la dens----
# # Parametres
# # Création d'une colonne ID
# data %>%
#   rowid_to_column("ID") %>%
#   filter(Site == "OVF") -> rf_data
#
# indices_rmf <- CreateSpacetimeFolds(rf_data,
#                                     spacevar = "Bloc",
#                                     k = 5,
#                                     seed = 2410)
#
# # Configuration de trainControl
# ctrl <- trainControl(method = "cv",
#                      index = indices_rmf$index,
#                      savePredictions = "final",
#                      summaryFunction = defaultSummary)
#
# #Hyperparameter tuning grid
# tuneGrid_rf <- expand.grid(.mtry = c(2, 5, 10, 15))
#
# set.seed(2410)
# num.trees_values = c(500, 2500, 5000, 10000, 15000, 20000)
# result <- data.frame()
#
# cl <- makeCluster(detectCores() - 1)
# registerDoParallel(cl)
#
# i <- 500
#
# for(i in num.trees_values){
#
#   print(i)
#
#   start.time = Sys.time()
#
#   set.seed(123)
#   train(x = rf_data[, selected_var],
#         y = rf_data$dens,
#         method = "rf",
#         trControl = ctrl,
#         tuneGrid = tuneGrid_rf,
#         metric = "RMSE",
#         ntree = i,
#         importance = TRUE) -> RF_1
#
#   varImp(RF_1)$importance %>%
#     as.data.frame() %>%
#     mutate(variable = row.names(imp_var_rf1)) -> VI_1
#
#   set.seed(321)
#
#   train(x = rf_data[, selected_var],
#         y = rf_data$dens,
#         method = "rf",
#         trControl = ctrl,
#         tuneGrid = tuneGrid_rf,
#         metric = "RMSE",
#         ntree = i,
#         importance = TRUE) -> RF_2
#
#   varImp(RF_2)$importance %>%
#     as.data.frame() %>%
#     mutate(variable = row.names(imp_var_rf1)) -> VI_2
#
#   end.time = Sys.time()
#
#   M <- merge(VI_1, VI_2, by="variable")
#   result = rbind(result, data.frame(number_of_trees = i,
#                                     variable_importance_stability = cor(M$Overall.x, M$Overall.y),
#                                     computation_time = (end.time - start.time)/2))
# }
#
# stopCluster(cl)
# cat("Training complete.\n")
#
# result
#
# optRF_result <- opt_importance(y = rf_data$dens, X = rf_data)
# summary(optRF_result) # utiliser le nombre recommande d'arbres pour la suite -> 5000
#
#
#
#

# # Selection des variables pour les modeles random forest par iteration PY ----
# p <- 5
# for(p in seq_len(nrow(params))){
#   param <- params[p,]
#   model_name <- paste0(param$site, "_", param$variable)
#
#   cat("\n--- compute random forest for model :", model_name, "---\n")
#
#   # Parametres
#   # Création d'une colonne ID
#   data %>%
#     rowid_to_column("ID") %>%
#     filter(Site == param$site) -> rf_data
#
#   indices_rmf <- CreateSpacetimeFolds(rf_data,
#                                       spacevar = "Bloc",
#                                       k = 5,
#                                       seed = 2410)
#
#   # Configuration de trainControl
#   ctrl <- trainControl(method = "cv",
#                        index = indices_rmf$index,
#                        savePredictions = "final",
#                        summaryFunction = defaultSummary)
#
#
#   # Old version PY
#   set.seed(2410)
#   remaining_vars <- selected_var
#   rf_models <- list()
#   model_rmse_perf <- data.frame()
#   model_r2_perf <- data.frame()
#   var_importance_list <- list()
#
#   cl <- makeCluster(detectCores() - 1)
#   registerDoParallel(cl)
#
#   for(i in 1:(length(selected_var) - 1)) {
#
#     rf_model <- train(x = rf_data[,remaining_vars],
#                       y = rf_data[,param$variable],
#                       method = "rf",
#                       trControl = ctrl,
#                       tuneGrid = tuneGrid_rf,
#                       metric = "RMSE",
#                       ntree = 7000,
#                       importance = TRUE)
#
#     rf_models[[i]] <- rf_model
#
#     imp <- varImp(rf_model)$importance %>%
#       as.data.frame() %>%
#       mutate(variable = row.names(.)) %>%
#       arrange(Overall)
#
#     var_importance_list[[i]] <- imp
#
#     least_important <- imp$variable[1]
#
#     model_rmse <- rf_model$results %>% mutate(iter = i)
#
#     model_r2 <- tail(rf_model$finalModel$rsq, 1) * 100 %>% as.data.frame()
#
#     model_r2$iter <- i
#
#     model_rmse_perf <- rbind(model_rmse_perf, model_rmse)
#
#     model_r2_perf <- rbind(model_r2_perf, model_r2)
#
#     cat("--- iteration", i,
#         "/ number of variable used :", length(remaining_vars),
#         "/ r2 :", model_r2[,1],
#         "/ rmse min :", min(model_rmse$RMSE),
#         "/ removing least important variable :", least_important, "---\n")
#
#      remaining_vars <- setdiff(remaining_vars, least_important)
#
#     if(length(remaining_vars) <= 2) break
#   }
#
#   stopCluster(cl)
#
#   cat("training complete\n")
#
#   var_imp <- enframe(var_importance_list, name = "iter", value = "values") %>%
#     unnest(values)
#
#   performance <- model_rmse_perf %>%
#     slice_min(RMSE) %>%
#     left_join(model_r2_perf, by = join_by(iter)) %>%
#     rename(R2 = '.') %>%
#     left_join(var_imp, by = join_by(iter)) %>%
#     arrange(-Overall) %>%
#     mutate(Model = model_name)
#
#   iter <- unique(performance$iter)
#
#   cat("best iteration :", iter,
#       "/ r2 :", unique(performance$R2),
#       "/ rmse :", unique(performance$RMSE),
#       "/ variables used :", paste0(performance$variable, collapse = " - "), "\n")
#
#   cat("save performance table\n")
#
#   write.csv(performance, paste0("./results/PY/performance_", model_name, ".csv"), row.names = FALSE)
#
#   final_model <- rf_models[[iter]]$finalModel
#
#   cat("save model\n")
#
#   saveRDS(final_model, paste0("./results/PY/model_", model_name, ".rds"))
#
# }
#
#
#
#

# Selection des variables pour les modeles random forest par iteration FL ----
s <- "RMF"
v <- "st"

for(s in sites){

  cat("build clusters and correlation matrix for site :", s, "\n")

  # Correlation entre les variables pour la selection
  data %>%
    filter(Site == s) %>%
    dplyr::select(all_of(selected_var)) %>%
    cor(use = "complete.obs") -> corr_matrix

  # # old selection method
  # selected_var[!selected_var %in% findCorrelation(corr_matrix, cutoff = threshold_corr, names = TRUE)] -> selected_var_filtered

  # Distance based on correlation
  dist_matrix <- as.dist(1 - abs(corr_matrix))

  # Hierarchical clustering
  hc <- hclust(dist_matrix, method = "average")

  clusters <- cutree(hc, h = 1 - threshold_corr)

  hc %>%
    cutree(h = 1 - threshold_corr) %>%
    tibble(var = names(.),
           cluster = .) %>%
    arrange(cluster) %>%
    rowid_to_column("id") -> cluster

  map_dfr(cluster$cluster %>% unique,
          function(x){

            cluster %>%
              filter(cluster == x) -> cluster_corr_temp

            tibble(xmin = min(cluster_corr_temp$id)-0.5,
                   xmax = max(cluster_corr_temp$id)+0.5,
                   ymin = min(cluster_corr_temp$id)-0.5,
                   ymax = max(cluster_corr_temp$id)+0.5)
          }) -> cluster_rectangles

  corr_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Variable1") %>%
    pivot_longer(-Variable1,
                 names_to = "Variable2",
                 values_to = "Correlation") %>%
    mutate(Variable1 = factor(Variable1, levels = cluster$var),
           Variable2 = factor(Variable2, levels = cluster$var)) -> corr_matrix_long

  # Create the heatmap plot with ggplot2
  ggplot(data = corr_matrix_long,
         mapping = aes(x = Variable1, y = Variable2, fill = Correlation)) +
    geom_tile(color = "white") +
    # new selection method
    geom_rect(mapping = aes(xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax),
              data = cluster_rectangles,
              inherit.aes = FALSE,
              color = "black", fill = NA, size = 1.5) +
    # # old selection method
    # geom_tile(data = corr_matrix_long %>% filter(Variable1 %in% selected_var_filtered & Variable2 %in% selected_var_filtered),
    #           mapping =  aes(Variable2, Variable1), color = "black", size = 2) + # Select filtered variables and display them with black rectangles
    geom_text(aes(label = round(Correlation, 2)),
              color = "black", size = 3) +
    geom_rect(size=1, fill=NA, colour="black",
              aes(xmin=0.5, xmax=0.5, ymin=0.5, ymax=0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limit = c(-1, 1), name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Reordered Correlation Matrix (Hierarchical Clustering)", x = "", y = "") -> correlation_matrix

  ggsave(paste0("./results/FL/", s, "_correlation_matrix.jpg"),
         plot = correlation_matrix,
         width = 20,
         height = 20,
         dpi = 300)

  # # old selection method
  # # Create the heatmap trimmed plot with ggplot2
  # ggplot(data = corr_matrix_long %>% filter(Variable1 %in% selected_var_filtered & Variable2 %in% selected_var_filtered),
  #        mapping = aes(x = Variable1, y = Variable2, fill = Correlation)) +
  #   geom_tile(color = "white") +
  #   geom_text(aes(label = round(Correlation, 2)),
  #             color = "black", size = 3) +
  #   geom_rect(size=1, fill=NA, colour="black",
  #             aes(xmin=0.5, xmax=0.5, ymin=0.5, ymax=0.5)) +
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
  #                        limit = c(-1, 1), name = "Correlation") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   labs(title = "Reordered Correlation Matrix (Hierarchical Clustering)", x = "", y = "")
  #
  # ggsave(paste0("./results/FL/", s, "_correlation_matrix_trim.jpg"),
  #        width=12,
  #        height=12,
  #        dpi=300)

  for(v in variables){

    model_name <- paste0(s, "_", v)

    cat("\ncompute random forest for model :", model_name, "\n")

    cat("model :", model_name, "\n\n",
        file = paste0("./results/FL/", model_name, "_log.txt"),
        append = TRUE)

    # Select site
    data %>%
      rowid_to_column("ID") %>%
      filter(Site == s) -> rf_data

    # Select the variable with highest pearson R per cluster
    map_vec(unique(cluster$cluster),
            function(x){

              cluster %>%
                filter(cluster == x) %>%
                pull(var) -> cluster_var

              cor(rf_data[,c(v, cluster_var)],
                  use = "complete.obs",
                  method = "pearson") %>%
                as_tibble() %>%
                mutate(name = colnames(.)) %>%
                dplyr::select(val = v,
                              name) %>%
                filter(name != v) %>%
                mutate(val = abs(val)) %>%
                arrange(desc(val)) %>%
                slice(1) %>%
                pull(name) -> name_selected

              return(name_selected)

            }) -> selected_var_filtered

    # Create the heatmap trimmed plot with ggplot2
    ggplot(data = corr_matrix_long %>% filter(Variable1 %in% selected_var_filtered & Variable2 %in% selected_var_filtered),
           mapping = aes(x = Variable1, y = Variable2, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(Correlation, 2)),
                color = "black", size = 3) +
      geom_rect(size=1, fill=NA, colour="black",
                aes(xmin=0.5, xmax=0.5, ymin=0.5, ymax=0.5)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           limit = c(-1, 1), name = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Reordered Correlation Matrix (Hierarchical Clustering)", x = "", y = "") -> correlation_matrix_trimmed

    ggsave(paste0("./results/FL/", model_name, "_correlation_matrix.jpg"),
           plot = correlation_matrix_trimmed,
           width=12,
           height=12,
           dpi=300)

    cat("all variables (cluster) :\n", paste0(cluster %>% mutate(label = paste0(var, " (", cluster, ")")) %>% pull(label), collapse = " - "), "\n\n",
        file = paste0("./results/FL/", model_name, "_log.txt"),
        append = TRUE)

    cat("filtered variables :\n", paste0(selected_var_filtered, collapse = " - "), "\n\n",
        file = paste0("./results/FL/", model_name, "_log.txt"),
        append = TRUE)

    cat("modelling method : random forest\nnumber of k-fold : 5\nmethod : repeatedcv\nrepeats : 10\nntree : 7000\nmetric : RMSE\nvariable importance method : permutation\nvariable importance wrapper prediction method : raw\nvariable importance metric : RMSE\nvariable importance nsim : 2\nseed : 1\n\n",
        file = paste0("./results/FL/", model_name, "_log.txt"),
        append = TRUE)

    # Creation des k-fold
    indices_rmf <- CreateSpacetimeFolds(rf_data,
                                        spacevar = "Bloc",
                                        k = 5,
                                        seed = 1)

    # Configuration de trainControl
    ctrl <- trainControl(method = "repeatedcv", # a valider
                         repeats = 10, # a valider
                         index = indices_rmf$index,
                         savePredictions = "final",
                         summaryFunction = defaultSummary)

    remaining_vars <- selected_var_filtered
    results_rf <- list()

    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    set.seed(1)

    while(length(remaining_vars) >= 2){

      set.seed(1)
      rf_model <- train(x = rf_data[,remaining_vars],
                        y = rf_data[,v],
                        method = "rf",
                        trControl = ctrl,
                        tuneGrid = tuneGrid_rf,
                        metric = "RMSE",
                        ntree = 7000, # ?????????????????????
                        importance = TRUE)

      set.seed(1)
      vip::vi(rf_model,
              method = "permute",
              train = rf_data[,c(remaining_vars, v)],
              target = v,
              metric = "rmse",
              #metric = "mae", # avoir si c'est mieux
              pred_wrapper = pred_wrapper,
              smaller_is_better = FALSE,
              nsim = 2) -> variable_importance

      # Version variable removed
      variable_importance %>%
        arrange(Importance) %>%
        slice(1) %>%
        pull(Variable) -> variable_removed

      mtry_selected <- unique(rf_model[["pred"]][["mtry"]])

      rf_model$results %>%
        filter(mtry == mtry_selected) -> rf_model_final_results

      list(model = rf_model,
           mtry = mtry_selected,
           rmse = rf_model_final_results$RMSE,
           mae = rf_model_final_results$MAE,
           rsquared = rf_model_final_results$Rsquared,
           r2 = tail(rf_model$finalModel$rsq, 1) * 100,
           variables_used = remaining_vars,
           variable_importance = variable_importance,
           variable_removed = variable_removed) %>%
        list() %>%
        append(results_rf) -> results_rf

      cat("number of variable used :", length(remaining_vars),
          "/ r2 :", tail(rf_model$finalModel$rsq, 1) * 100,
          "/ rmse :", rf_model_final_results$RMSE,
          "/ mae :", rf_model_final_results$MAE,
          "/ variable removed :", variable_removed, "\n")
      cat("number of variable used :", length(remaining_vars),
          "/ r2 :", tail(rf_model$finalModel$rsq, 1) * 100,
          "/ rmse :", rf_model_final_results$RMSE,
          "/ mae :", rf_model_final_results$MAE,
          "/ variable removed :", variable_removed, "\n",
          file = paste0("./results/FL/", model_name, "_log.txt"),
          append = TRUE)

      remaining_vars <- setdiff(remaining_vars, variable_removed)

    }

    stopCluster(cl)

    results_rf %>%
      map_dfr(function(x){

        tibble(model = list(x$model),
               r2 = x$r2,
               rmse = x$rmse,
               mae = x$mae,
               n_variables_used = length(x$variables_used),
               variables_used = paste0(unlist(x$variables_used), collapse = " - "))

      }) -> results_rf_simplify

    scale_factor <- max(results_rf_simplify$r2) / max(results_rf_simplify$rmse)

    ggplot(data = results_rf_simplify,
           aes(x = n_variables_used)) +
      geom_line(aes(y = r2), color = "green2") +
      geom_line(aes(y = rmse * scale_factor), color = "red") +
      scale_x_reverse() +
      scale_y_continuous(name = "R2 (green)",
                         sec.axis = sec_axis(~ . / scale_factor, name = "RMSE (red)")) +
      labs(x = "N variable") -> plot_performance

    ggsave(paste0("./results/FL/", model_name, "_plot_performance.jpg"),
           plot = plot_performance,
           width = 10,
           height = 5,
           dpi = 300,
           units = "in")

    results_rf_simplify %>%
      #arrange(desc(r2)) %>%
      arrange(rmse) %>%
      slice(1) -> results_rf_best_model

    cat("best model params :",
        " r2 :", results_rf_best_model$r2,
        "/ rmse :", results_rf_best_model$rmse,
        "/ mae :", results_rf_best_model$mae,
        "/ variables used :", results_rf_best_model$variables_used, "\n")
    cat("\nbest model params :",
        " r2 :", results_rf_best_model$r2,
        "/ rmse :", results_rf_best_model$rmse,
        "/ mae :", results_rf_best_model$mae,
        "/ variables used :", results_rf_best_model$variables_used, "\n\n",
        file = paste0("./results/FL/", model_name, "_log.txt"),
        append = TRUE)

    cat("save performance table\n")

    results_rf_best_model %>%
      dplyr::select(-model) %>%
      write.csv(paste0("./results/FL/", model_name, "_performance.csv"), row.names = FALSE)

    cat("save model\n")

    results_rf_best_model %>%
      pull(model) %>%
      {.[[1]]$finalModel} %>%
      saveRDS(paste0("./results/FL/", model_name, "_model.rds"))

  }

}





# Validation des modèles ----
for(p in seq_len(nrow(params))){

  param <- params[p,]
  model_name <- paste0(param$site, "_", param$variable)

  readRDS(paste0("./results/FL/model_", model_name, ".rds")) -> final_model_FL

  readRDS(paste0("D:/00_Ontario_eFRI/random_forest/dossier_remise_fin_contrat_PY/modeles/final_mod_", param$variable, "_", param$site, ".rds")) -> final_model_PY

  cat("modele :", model_name, "\n",
      "/ var FL :", paste0(final_model_FL[["xNames"]], collapse = " - "), "\n",
      "/ var PY :", paste0(final_model_PY[["xNames"]], collapse = " - "), "\n")
}





# Application des modeles ----
modeles <- list.files("D:/dossier_remise_fin_contrat_PY/ontario/modeles", full.names = T)

dfa <- sf::st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp")

# Apply models
# m <- modeles[1]
for(m in modeles){

  modele_info <- tools::file_path_sans_ext(basename(m))
  modele_info <- stringr::str_remove(modele_info, "final_mod_")
  modele_info <- stringr::str_split(modele_info, "_")[[1]]
  dendrometric <- modele_info[[1]]
  forest <- toupper(modele_info[[2]])

  sectors <- dfa[dfa$forest == forest, ]$sector

  modele <- readRDS(m)
  metric_names <- modele[["xNames"]]
  metric_names <- gsub("^zq", "z_p", metric_names)
  metric_names <- gsub("\\.", "-", metric_names)
  metric_names <- paste0(metric_names, ".tif")
  metric_names <- gsub("dem.tif", "dem_buffered.tif", metric_names)
  metric_names <- gsub("slope.tif", "slope_buffered.tif", metric_names)
  metric_names <- gsub("0-15", "0.15", metric_names)
  metric_names <- gsub("pzabovemean", "z_abovemean", metric_names)

  # s <- sectors[1]
  for(s in sectors){

    cat(dendrometric, forest, s, "\n")

    files <- list.files(paste0("D:/dossier_remise_fin_contrat_PY/ontario/metriques/", s, "/metrics"),
                        pattern = "\\.tif$",
                        full.names = T)

    files <- files[stringr::str_detect(files, paste(metric_names, collapse = "|"))]

    if(length(files) == length(metric_names)){

      metrics <- purrr::map(files, terra::rast)
      metrics <- purrr::map(metrics, function(x) {terra::resample(x, metrics[[1]])})
      metrics <- terra::rast(metrics)

      metric_names_raster <- names(metrics)
      metric_names_raster <- gsub("z_p", "zq", metric_names_raster)
      metric_names_raster <- gsub("-", ".", metric_names_raster)
      metric_names_raster <- gsub("dem_buffered", "dem", metric_names_raster)
      metric_names_raster <- gsub("slope_buffered", "slope", metric_names_raster)
      metric_names_raster <- gsub("z_abovemean", "pzabovemean", metric_names_raster)
      names(metrics) <- metric_names_raster

      prediction <- terra::predict(metricsABA[c("z_below_0")], modele, na.rm = TRUE)
      prediction <- terra::predict(metrics[c("z_below_0")], modele, na.rm = TRUE)

      terra::writeRaster(prediction, paste0("D:/dossier_remise_fin_contrat_PY/ontario/predictions/2026_05_14/", dendrometric, "_", forest, "_", s, ".tif"))

    } else {

      cat("Not all metrics are available\n")

    }

  }

  # Merge models
  predictions_files <- list.files("D:/dossier_remise_fin_contrat_PY/ontario/predictions/2026_05_14", pattern = "\\.tif$", full.names = TRUE)
  predictions_files <- stringr::str_subset(predictions_files, dendrometric)
  predictions_files <- stringr::str_subset(predictions_files, forest)
  predictions <- purrr::map(predictions_files, terra::rast)
  predictions <- sprc(predictions)
  predictions <- terra::merge(predictions)
  terra::writeRaster(predictions, paste0("D:/dossier_remise_fin_contrat_PY/ontario/predictions/2026_05_14/", dendrometric, "_", forest, ".tif"))

}
