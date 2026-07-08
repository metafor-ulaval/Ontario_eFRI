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
library(tidyterra)
library(randomForest)
library(plotly)
library(CAST)
library(caret)
library(doParallel)
library(optRF)
library(pheatmap)
library(sf)
library(terra)
library(vip)

filter_lidar_parameters <- function(data,
                                    radius = "1128",
                                    zmin = "NA"){
  data %>%
    dplyr::select(Name, contains(paste0("_radius_", radius, "_cm"))) %>%
    dplyr::rename_with(.fn = ~ gsub(paste0("_radius_", radius, "_cm"), "", .)) %>%
    dplyr::select(Name, contains(paste0("_zmin_", zmin, "_cm"))) %>%
    dplyr::rename_with(.fn = ~ gsub(paste0("_zmin_", zmin, "_cm"), "", .)) -> data_filtered

  return(data_filtered)
}

filter_ufc_parameters <- function(data,
                                  radius = "1128"){
  data %>%
    dplyr::select(Name, contains(paste0("_radius_", radius))) %>%
    dplyr::rename_with(.fn = ~ gsub(paste0("_radius_", radius), "", .)) -> data_filtered

  return(data_filtered)
}

correct_placettes_names <- function(data){

  data %>%
    dplyr::mutate(Name = gsub("_centre", "", Name)) %>% # enleve l'appelation centre
    dplyr::mutate(Name = gsub("X", "", Name)) %>% # enleve les X dans les noms de parcelles qui indiquent que la parcelle a due etre deplacee
    dplyr::filter(!grepl('FC', Name)) %>% # enleve la parcelle FC
    dplyr::mutate(Name = case_when(Name == 'S10_52' ~ 'S10_5', .default = Name)) -> data_corrected # corrige un nom de parcelle

  return(data_corrected)
}

rename_pz <- function(data){

  data %>%
    dplyr::rename("pz_0-0.15" = "pz_0.0.15",
                  "pz_0.15-2" = "pz_0.15.2",
                  "pz_2-5" = "pz_2.5",
                  "pz_5-10" = "pz_5.10",
                  "pz_10-20" = "pz_10.20",
                  "pz_20-30" = "pz_20.30") -> data_renamed

  return(data_renamed)

}

# Les plus complexes qui ont ete mise de cotes pour faciliter le calcul des metriques
selected_var <- c("zcv",
                  # "zskew", # Retire a cause du pattern carre
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
                  # "zq95", # Retire a cause du pattern carre
                  # "zq99", # Retire a cause du pattern carre
                  "pzabovemean",
                  "pzabove2",
                  "pzabove5",
                  # "zpcum1", # Retire a cause du pattern carre
                  # "zpcum2", # Retire a cause du pattern carre
                  # "zpcum3", # Retire a cause du pattern carre
                  # "zpcum4", # Retire a cause du pattern carre
                  # "zpcum5", # Retire a cause du pattern carre
                  # "zpcum6", # Retire a cause du pattern carre
                  # "zpcum7", # Retire a cause du pattern carre
                  # "zpcum8", # Retire a cause du pattern carre
                  # "zpcum9", # Retire a cause du pattern carre
                  # "pz_below_0", # Retire a cause du pattern carre
                  "pz_0-0.15",
                  "pz_0.15-2",
                  "pz_2-5",
                  "pz_5-10",
                  "pz_10-20",
                  # "pz_20-30", # Retire a cause du pattern carre
                  # "pz_above_30", # pz_above_30 elimine car ne contient que des 0 pour RMF
                  "ziqr",
                  # "zMADmean", # Retire a cause du pattern carre
                  "zMADmedian",
                  "CRR",
                  # "rumple", # pas genere en buffered
                  # "zentropy", # force d'enlever zentropy car valeur manquante pour une des parcelles improductives de OVF
                  # "ufc_circle_1m_to_2m", # pas genere en buffered
                  # "ufc_circle_2m_to_3m", # pas genere en buffered
                  # "ufc_circle_3m_to_4m", # pas genere en buffered
                  # "ufc_circle_1m_to_4m", # pas genere en buffered
                  # "sagawi",
                  # "twi",
                  "slope",
                  "eastness",
                  "northness",
                  # "tpi",
                  "dem")

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
readRDS("./base_data/placettes_data_topo.rds") %>%
  dplyr::select(-sector) %>%
  correct_placettes_names() %>%
  rename_with(~ sub("_mean$", "", .)) -> topo


# Metriques lidar ufc
readRDS("./base_data/placettes_data_improductives.rds") %>%
  st_drop_geometry() %>%
  filter_ufc_parameters() %>%
  dplyr::select(Name, contains("ufc")) -> imp_ufc

readRDS("./base_data/placettes_data.rds") %>%
  filter_ufc_parameters() %>%
  dplyr::select(Name, contains("ufc")) %>%
  correct_placettes_names() %>%
  bind_rows(imp_ufc) -> ufc


# Metriques lidar set3
readRDS("./base_data/placettes_data_improductives.rds") %>%
  st_drop_geometry() %>%
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

rm(topo,
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
  mutate(vmb_ha = vmb * 25) %>%
  full_join(data, by = "Name") %>%
  filter(Name != "S7_4") %>% # enlève une ligne avec un nom de parcelle probablement errone (sans match dans les metriques)
  filter(Name != "S5_3") %>% # enlève une parcelle ou les donnes lidar sont erronees
  mutate(Site = case_when(Bloc %in% paste0("S", 1:5) ~ "RMF",
                          Bloc %in% paste0("S", 6:10) ~ "OVF")) -> data





# Selection des variables pour les modeles random forest par iteration ----
# s <- "RMF"
# v <- "dens"

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
              color = "black", fill = NA, linewidth = 1.5) +
    # # old selection method
    # geom_tile(data = corr_matrix_long %>% filter(Variable1 %in% selected_var_filtered & Variable2 %in% selected_var_filtered),
    #           mapping =  aes(Variable2, Variable1), color = "black", size = 2) + # Select filtered variables and display them with black rectangles
    geom_text(aes(label = round(Correlation, 2)),
              color = "black", size = 3) +
    geom_rect(linewidth=1, fill=NA, colour="black",
              aes(xmin=0.5, xmax=0.5, ymin=0.5, ymax=0.5)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limit = c(-1, 1), name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Reordered Correlation Matrix (Hierarchical Clustering)", x = "", y = "") -> correlation_matrix

  ggsave(paste0("./results/", s, "_correlation_matrix.jpg"),
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
  # ggsave(paste0("./results/", s, "_correlation_matrix_trim.jpg"),
  #        width=12,
  #        height=12,
  #        dpi=300)

  for(v in variables){

    model_name <- paste0(s, "_", v)

    cat("\ncompute random forest for model :", model_name, "\n")

    cat("model :", model_name, "\n\n",
        file = paste0("./results/", model_name, "_log.txt"),
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
                dplyr::select(val = all_of(v),
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
      geom_rect(linewidth=1, fill=NA, colour="black",
                aes(xmin=0.5, xmax=0.5, ymin=0.5, ymax=0.5)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           limit = c(-1, 1), name = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Reordered Correlation Matrix (Hierarchical Clustering)", x = "", y = "") -> correlation_matrix_trimmed

    ggsave(paste0("./results/", model_name, "_correlation_matrix.jpg"),
           plot = correlation_matrix_trimmed,
           width=12,
           height=12,
           dpi=300)

    # Creation des k-fold
    indices_rmf <- CreateSpacetimeFolds(rf_data,
                                        spacevar = "Bloc",
                                        k = 5,
                                        seed = 1)

    # Configuration de trainControl
    ctrl <- trainControl(method = "repeatedcv",
                         repeats = 10,
                         index = indices_rmf$index,
                         savePredictions = "final",
                         summaryFunction = defaultSummary)

    remaining_vars <- selected_var_filtered
    results_rf <- list()

    # Find the optimal number of trees
    set.seed(1)
    opt_importance(y = rf_data[,v],
                   X = rf_data[,remaining_vars],
                   num.trees_values = c(250, 500, 750, 1000, 2000, 5000, 10000)) -> opt_ntrees_results
    opt_ntrees_results$recommendation -> opt_ntrees

    cat("all variables (cluster) :\n", paste0(cluster %>% mutate(label = paste0(var, " (", cluster, ")")) %>% pull(label), collapse = " - "), "\n\n",
        file = paste0("./results/", model_name, "_log.txt"),
        append = TRUE)

    cat("correlation threshold for clustering and filtering based on pearson correlation :", threshold_corr, "\n\n",
        file = paste0("./results/", model_name, "_log.txt"),
        append = TRUE)

    cat("filtered variables :\n", paste0(selected_var_filtered, collapse = " - "), "\n\n",
        file = paste0("./results/", model_name, "_log.txt"),
        append = TRUE)

    cat(paste0("modelling method : random forest\nnumber of k-fold : 5\nmethod : repeatedcv\nrepeats : 10\nntree : ", opt_ntrees, "\nmetric : RMSE\nvariable importance method : permutation\nvariable importance wrapper prediction method : raw\nvariable importance metric : RMSE\nvariable importance nsim : 10\nseed : 1\n\n"),
        file = paste0("./results/", model_name, "_log.txt"),
        append = TRUE)

    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)

    set.seed(1)

    while(length(remaining_vars) >= 3){

      set.seed(1)
      rf_model <- train(x = rf_data[,remaining_vars],
                        y = rf_data[,v],
                        method = "rf",
                        trControl = ctrl,
                        tuneGrid = tuneGrid_rf,
                        metric = "RMSE",
                        ntree = opt_ntrees,
                        importance = TRUE)

      set.seed(1)
      vip::vi(rf_model,
              method = "permute",
              train = rf_data[,c(remaining_vars, v)],
              target = v,
              metric = "rmse",
              #metric = "mae", # a voir si c'est mieux
              pred_wrapper = pred_wrapper,
              smaller_is_better = FALSE,
              nsim = 10) -> variable_importance

      # Version variable removed
      variable_importance %>%
        arrange(Importance) %>%
        slice(1) %>%
        pull(Variable) -> variable_removed

      mtry_selected <- unique(rf_model[["pred"]][["mtry"]])

      rf_model$results %>%
        filter(mtry == mtry_selected) -> rf_model_final_results

      rf_data %>%
        mutate(pred = predict(rf_model, .),
               error = pred - !!sym(v)) %>%
        pull(error) -> errors

      list(model = rf_model,
           mtry = mtry_selected,
           rmse = rf_model_final_results$RMSE,
           mae = rf_model_final_results$MAE,
           rsquared = rf_model_final_results$Rsquared,
           r2 = tail(rf_model$finalModel$rsq, 1) * 100,
           bias = mean(errors),
           variables_used = remaining_vars,
           variable_importance = variable_importance,
           variable_removed = variable_removed) %>%
        list() %>%
        append(results_rf) -> results_rf

      cat("number of variable used :", length(remaining_vars),
          "/ r2 :", tail(rf_model$finalModel$rsq, 1) * 100,
          "/ rmse :", rf_model_final_results$RMSE,
          "/ mae :", rf_model_final_results$MAE,
          "/ bias :", mean(errors),
          "/ variable removed :", variable_removed, "\n")
      cat("number of variable used :", length(remaining_vars),
          "/ r2 :", tail(rf_model$finalModel$rsq, 1) * 100,
          "/ rmse :", rf_model_final_results$RMSE,
          "/ mae :", rf_model_final_results$MAE,
          "/ bias :", mean(errors),
          "/ variable removed :", variable_removed, "\n",
          file = paste0("./results/", model_name, "_log.txt"),
          append = TRUE)

      remaining_vars <- setdiff(remaining_vars, variable_removed)

    }

    stopCluster(cl)

    results_rf %>%
      saveRDS(paste0("./results/", model_name, "_results_all.rds"))

    results_rf %>%
      map_dfr(function(x){

        tibble(model = list(x$model),
               r2 = x$r2,
               rmse = x$rmse,
               mae = x$mae,
               bias = x$bias,
               n_variables_used = length(x$variables_used),
               variables_used = paste0(unlist(x$variables_used), collapse = " - "))

      }) -> results_rf_simplify

    results_rf_simplify %>%
      saveRDS(paste0("./results/", model_name, "_results_simplify.rds"))

    scale_factor_rmse <- max(results_rf_simplify$r2) / max(results_rf_simplify$rmse)

    ggplot(data = results_rf_simplify,
           aes(x = n_variables_used)) +
      geom_line(aes(y = r2), color = "green2") +
      geom_line(aes(y = rmse * scale_factor_rmse), color = "red") +
      scale_x_reverse() +
      scale_y_continuous(name = "R2 (green)",
                         sec.axis = sec_axis(~ . / scale_factor_rmse, name = "RMSE (red)")) +
      labs(x = "N variable") -> plot_performance_rmse

    ggsave(paste0("./results/", model_name, "_plot_performance_rmse.jpg"),
           plot = plot_performance_rmse,
           width = 10,
           height = 5,
           dpi = 300,
           units = "in")

    scale_factor_mae <- max(results_rf_simplify$r2) / max(results_rf_simplify$mae)

    ggplot(data = results_rf_simplify,
           aes(x = n_variables_used)) +
      geom_line(aes(y = r2), color = "green2") +
      geom_line(aes(y = mae * scale_factor_mae), color = "blue2") +
      scale_x_reverse() +
      scale_y_continuous(name = "R2 (green)",
                         sec.axis = sec_axis(~ . / scale_factor_mae, name = "MAE (blue)")) +
      labs(x = "N variable") -> plot_performance_mae

    ggsave(paste0("./results/", model_name, "_plot_performance_mae.jpg"),
           plot = plot_performance_mae,
           width = 10,
           height = 5,
           dpi = 300,
           units = "in")

    results_rf_simplify %>%
      #arrange(desc(r2)) %>% # a voir si c'est mieux
      arrange(rmse) %>%
      slice(1) -> results_rf_best_model

    cat("best model params :",
        " r2 :", results_rf_best_model$r2,
        "/ rmse :", results_rf_best_model$rmse,
        "/ mae :", results_rf_best_model$mae,
        "/ bias :", results_rf_best_model$bias,
        "/ variables used :", results_rf_best_model$variables_used, "\n")
    cat("\nbest model params :",
        " r2 :", results_rf_best_model$r2,
        "/ rmse :", results_rf_best_model$rmse,
        "/ mae :", results_rf_best_model$mae,
        "/ bias :", results_rf_best_model$bias,
        "/ variables used :", results_rf_best_model$variables_used, "\n\n",
        file = paste0("./results/", model_name, "_log.txt"),
        append = TRUE)

    cat("save performance table\n")

    results_rf_best_model %>%
      dplyr::select(-model) %>%
      write.csv(paste0("./results/", model_name, "_performance.csv"), row.names = FALSE)

    cat("save model\n")

    results_rf_best_model %>%
      pull(model) %>%
      {.[[1]]$finalModel} %>%
      saveRDS(paste0("./results/", model_name, "_model.rds"))

  }

}





# Creation d'une table finale de resultats ----
list.files("./results", pattern = "results_simplify.rds", full.names = T) %>%
  map_dfr(function(x){

    modele_info <- tools::file_path_sans_ext(basename(x)) %>% str_remove("_results_simplify")
    dendrometric <- substr(modele_info, 5, nchar(modele_info))
    forest <- substr(modele_info, 1, 3)

    readRDS(x) %>%
      dplyr::select(-model) %>%
      arrange(rmse) %>%
      slice(1) %>%
      mutate(forest = forest,
             dendrometric = dendrometric)
  }) %>%
  write.csv("./results/results_final.csv")






# Application des modeles ----
modeles <- list.files("./results/", pattern = "model", full.names = T)

dfa <- sf::st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp")

# Apply models
# m <- modeles[1]

for(m in modeles){

  modele <- readRDS(m)

  modele_info <- tools::file_path_sans_ext(basename(m)) %>% str_remove("_model")
  dendrometric <- substr(modele_info, 5, nchar(modele_info))
  forest <- substr(modele_info, 1, 3)
  sectors <- dfa[dfa$forest == forest, ]$sector
  metric_names <- modele[["xNames"]]

  # s <- sectors[1]

  for(s in sectors){

    cat(dendrometric, forest, s, "\n")

    paste0("D:/00_Ontario_eFRI/random_forest/metriques/", s, "/metrics/", metric_names, ".tif") %>%
      map(rast) %>%
      map(resample, .[[1]]) %>%
      rast() -> metrics

    names(metrics) <- metric_names

    metrics %>%
      terra::predict(modele, na.rm = TRUE) -> prediction

    names(prediction) <- dendrometric

    ggplot() +
      geom_spatraster(data = prediction) +
      labs(title = paste0(forest, " / ", s),
           fill = dendrometric) +
      tidyterra::scale_fill_whitebox_c(palette = "bl_yl_rd") -> prediction_plot

    ggsave(paste0("./results/", forest, "_", dendrometric, "_", s, "_prediction_plot.jpg"),
           plot = prediction_plot,
           width=12,
           height=12,
           dpi=300)

    prediction %>%
      terra::writeRaster(paste0("./results/", forest, "_", dendrometric, "_", s, "_prediction.tif"))

  }

  # Merge models
  list.files("./results/", pattern = "\\.tif$", full.names = TRUE) %>%
    str_subset(paste0(forest, "_", dendrometric)) %>%
    map(rast) %>%
    sprc() %>%
    terra::merge() %>%
    writeRaster(paste0("./results/", forest, "_", dendrometric, "_prediction.tif"))

}
