# Questions
# Pourquoi on fait une corrélation cutoff de 0.9 pour connaitre les variables a éléminer, mais au final on garde tout ?
# Pourquoi on fait deux RF pour le nombre optimal d'arbres ?
# Pourquoi on fait le nombre optimal d'arbres alors qu'il existe une fonction intégrée qui le fait plus rapidement ? (opt_importance())
# Pourquoi 7000 arbres alors que la fonction et l'analyse de sensibilité en propose 5000 ?
# Pourquoi ne pas mettre le set.seed juste avant le train ou bien dans le ctrl du train ?
# Pourquoi les arbres ne sont pas faits dans ranger ?
# Pourquoi faire un k-cross-fold ? -> tout le monde fait ça
# Est-ce que la permutation a été envisagée pour retirer itérativement les variables du modèles ? La simple importance ne fonctionne pas lorsqu'il existe plusieurs variable corrélés entre elles

# PCA avec corrélation avec la variable réponse (prcomp)

# Johanne White 2017 - Propose une méthode?




# Librarie, functions and parameters ----
library(tidyverse)
library(randomForest)
library(CAST)
library(caret)
library(doParallel)
library(optRF)
library(terra)

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

sites <- c("RMF", "OVF")

variables <- c("vmb_ha", "st", "dhpq", "dens")

params <- expand.grid(site = sites,
                      variable = variables,
                      stringsAsFactors = FALSE)




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

# Metriques lidar
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
  full_join(topo, by = "Name") -> data

rm(imp_topo,
   topo,
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





# Correlation entre les variables ----
# s <- "RMF"
for(s in sites){
  print(s)
  data %>%
    filter(Site == s) %>%
    dplyr::select(all_of(selected_var)) %>%
    cor() %>%
    findCorrelation(cutoff = 0.9, name = TRUE) %>%
    print()
}





# Test nombre optimal d'arbres pour ovf et la dens----
# Parametres
# Création d'une colonne ID
data %>%
  rowid_to_column("ID") %>%
  filter(Site == "OVF") -> rf_data

indices_rmf <- CreateSpacetimeFolds(rf_data,
                                    spacevar = "Bloc",
                                    k = 5,
                                    seed = 2410)

# Configuration de trainControl
ctrl <- trainControl(method = "cv",
                     index = indices_rmf$index,
                     savePredictions = "final",
                     summaryFunction = defaultSummary)

#Hyperparameter tuning grid
tuneGrid_rf <- expand.grid(.mtry = c(2, 5, 10, 15))

set.seed(2410)
num.trees_values = c(500, 2500, 5000, 10000, 15000, 20000)
result <- data.frame()

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

i <- 500

for(i in num.trees_values){

  print(i)

  start.time = Sys.time()

  set.seed(123)
  train(x = rf_data[, selected_var],
        y = rf_data$dens,
        method = "rf",
        trControl = ctrl,
        tuneGrid = tuneGrid_rf,
        metric = "RMSE",
        ntree = i,
        importance = TRUE) -> RF_1

  varImp(RF_1)$importance %>%
    as.data.frame() %>%
    mutate(variable = row.names(imp_var_rf1)) -> VI_1

  set.seed(321)

  train(x = rf_data[, selected_var],
        y = rf_data$dens,
        method = "rf",
        trControl = ctrl,
        tuneGrid = tuneGrid_rf,
        metric = "RMSE",
        ntree = i,
        importance = TRUE) -> RF_2

  varImp(RF_2)$importance %>%
    as.data.frame() %>%
    mutate(variable = row.names(imp_var_rf1)) -> VI_2

  end.time = Sys.time()

  M <- merge(VI_1, VI_2, by="variable")
  result = rbind(result, data.frame(number_of_trees = i,
                                    variable_importance_stability = cor(M$Overall.x, M$Overall.y),
                                    computation_time = (end.time - start.time)/2))
}

stopCluster(cl)
cat("Training complete.\n")

result

optRF_result <- opt_importance(y = rf_data$dens, X = rf_data)
summary(optRF_result) # utiliser le nombre recommande d'arbres pour la suite -> 5000





# Selection des variables pour les modeles random forest par iteration ----
for(p in seq_len(nrow(params))){
  param <- params[p,]
  model_name <- paste0(param$site, "_", param$variable)

  cat("\n--- compute random forest for model :", model_name, "---\n")

  # Parametres
  # Création d'une colonne ID
  data %>%
    rowid_to_column("ID") %>%
    filter(Site == param$site) -> rf_data

  indices_rmf <- CreateSpacetimeFolds(rf_data,
                                      spacevar = "Bloc",
                                      k = 5,
                                      seed = 2410)

  # Configuration de trainControl
  ctrl <- trainControl(method = "cv",
                       index = indices_rmf$index,
                       savePredictions = "final",
                       summaryFunction = defaultSummary)

  # À voir où mettre ces élément pour réduire la taille
  set.seed(2410)
  remaining_vars <- selected_var
  rf_models <- list()
  model_rmse_perf <- data.frame()
  model_r2_perf <- data.frame()
  var_importance_list <- list()
  # -------

  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)

  for(i in 1:(length(selected_var) - 1)) { # Remplacer par un while()

    rf_model <- train(x = rf_data[,remaining_vars],
                      y = rf_data[,param$variable],
                      method = "rf",
                      trControl = ctrl,
                      tuneGrid = tuneGrid_rf,
                      metric = "RMSE",
                      ntree = 7000,
                      importance = TRUE)

    rf_models[[i]] <- rf_model

    imp <- varImp(rf_model)$importance %>%
      as.data.frame() %>%
      mutate(variable = row.names(.)) %>%
      arrange(Overall)

    var_importance_list[[i]] <- imp

    least_important <- imp$variable[1]

    remaining_vars <- setdiff(remaining_vars, least_important)

    model_rmse <- rf_model$results %>% mutate(iter = i)

    model_r2 <- tail(rf_model$finalModel$rsq, 1) * 100 %>% as.data.frame()

    model_r2$iter <- i

    model_rmse_perf <- rbind(model_rmse_perf, model_rmse)

    model_r2_perf <- rbind(model_r2_perf, model_r2)

    cat("--- iteration", i,
        "/ number of variable used :", length(remaining_vars),
        "/ r2 :", model_r2[,1],
        "/ rmse min :", min(model_rmse$RMSE),
        "/ removing least important variable :", least_important, "---\n")

    if(length(remaining_vars) <= 2) break
  }

  stopCluster(cl)

  cat("training complete\n")

  var_imp <- enframe(var_importance_list, name = "iter", value = "values") %>%
    unnest(values)

  performance <- model_rmse_perf %>%
    slice_min(RMSE) %>%
    left_join(model_r2_perf, by = join_by(iter)) %>%
    rename(R2 = '.') %>%
    left_join(var_imp, by = join_by(iter)) %>%
    arrange(-Overall) %>%
    mutate(Model = model_name)

  iter <- unique(performance$iter)

  cat("best iteration :", iter,
      "/ r2 :", unique(performance$R2),
      "/ rmse :", unique(performance$RMSE),
      "/ variables used :", paste0(performance$variable, collapse = " - "), "\n")

  cat("save performance table\n")

  write.csv(performance, paste0("./results/performance_", model_name, ".csv"), row.names = FALSE)

  final_model <- rf_models[[iter]]$finalModel

  cat("save model\n")

  saveRDS(final_model, paste0("./results/model_", model_name, ".rds"))

}





# Validation des modèles ----
for(p in seq_len(nrow(params))){

  param <- params[p,]
  model_name <- paste0(param$site, "_", param$variable)

  readRDS(paste0("./results/model_", model_name, ".rds")) -> final_model_FL

  readRDS(paste0("D:/dossier_remise_fin_contrat_PY/ontario/modeles/final_mod_", param$variable, "_", param$site, ".rds")) -> final_model_PY

  cat("modele :", model_name,
      "/ var FL :", paste0(final_model_FL[["xNames"]], collapse = " - "),
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
