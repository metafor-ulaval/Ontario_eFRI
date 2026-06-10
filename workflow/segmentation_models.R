# Parameters ----
library(tidyverse)
library(magrittr)
library(terra)
library(tidyterra)
library(smoothr)
library(exactextractr)
library(rmapshaper)
library(sf)
library(pfif)
# remotes::install_github("metafor-ulaval/eFRItools")
library(eFRItools)
library(knitr)
library(basemapR)





# Parameters ----
otb_dir <- "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin"
Sys.setenv(OTB_MAX_RAM_HINT = "65536") # 64 GO
Sys.setenv(OTB_MEMORY_AVAILABLE = "65536") # 64 GO
Sys.setenv(GDAL_CACHEMAX = "65536") # 64 GO





# Models ----
# # First ones (keep them?)
# models <- list(model1 = c("z_p95", "z_cv", "z_above2"),
#                model2 = c("z_p95", "z_cv", "z_above2", "slope"),
#                model3 = c("z_p95", "z_cv", "z_above2", "fractional_cover_05_2"),
#                model4 = c("z_p95", "z_cv", "z_above2", "sagawi", "fractional_cover_05_2"),
#                model5 = c("z_p95", "z_cv", "z_above2", "NDMI"),
#                model6 = c("z_p95", "slope", "fractional_cover_05_2", "NDVI"),
#                model7 = c("z_p80", "sagawi", "z_skew", "DVI"),
#                model8 = c("z_p80", "slope", "fractional_cover_05_2", "NDMI"),
#                model9 = c("z_p80", "z_cv", "z_above2", "fractional_cover_05_2"),
#                model10 = c("z_p80", "z_cv", "sagawi", "fractional_cover_05_2"),
#                model11 = c("z_p80", "DVI", "sagawi", "z_above2"),
#                model12 = c("z_p80", "sagawi", "z_above2", "B4"))
#
#
# # First ones with optimised parameters (keep them?)
# # RMF
# wd <- "D:/00_Ontario_eFRI/RMF"
# epsg <- "EPSG:2958"
#
# models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
#                model_B = c("z_p95", "z_above2", "z_skew"),
#                model_C = c("z_p95", "z_above2", "fractional_cover_2_5"),
#                model_D = c("z_p95", "z_above2", "sagawi"),
#                model_E = c("z_p95", "z_above2", "slope"),
#                model_F = c("z_p95", "z_above2", "B6"),
#                model_G = c("B6", "B8", "B11"),
#                model_H = c("B2", "B3", "B4"),
#                model_I = c("V_ha", "dens", "qmdbh"))
#
# thresh <- 50
# spec <- 0.9
# spat <- 0.1
#
# # OVF
# wd <- "D:/00_Ontario_eFRI/OVF"
# epsg <- "EPSG:2959"
#
# models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
#                model_B = c("z_p95", "z_above2", "z_skew"),
#                model_D = c("z_p95", "z_above2", "sagawi"),
#                model_E = c("z_p95", "z_above2", "slope"),
#                model_F = c("z_p95", "z_above2", "B6"),
#                model_G = c("B6", "B8", "B11"),
#                model_H = c("B2", "B3", "B4"))
#
# thresh <- 60
# spec <- 0.5
# spat <- 1
#
# Models with optimised parameters
readRDS("D:/grm_parametrization/results/best_parameters_all.rds") %>%
  filter(type == "by_forest_model") %>%
  mutate(epsg = case_when(forest == "RMF" ~ "EPSG:2958",
                          forest == "OVF" ~ "EPSG:2959"))-> best_parameters_all

models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
               model_B = c("z_p95", "z_above2", "z_skew"),
               model_D = c("z_p95", "z_above2", "sagawi"),
               model_F = c("z_p95", "z_above2", "B6"),
               model_G = c("B6", "B8", "B11"))



# Computing
f <- "RMF"
m <- "model_B"

for(f in best_parameters_all$forest %>% unique()){
  cat(paste0(f, "\n"))
  wd <- paste0("D:/00_Ontario_eFRI/", f)

  # Read metrics
  list.files(paste0(wd, "/metrics"), full.names = TRUE) %>%
    read_metrics() -> metrics_infos

  best_parameters_all %>%
    filter(forest == f) %>%
    pull(epsg) %>%
    unique() -> epsg

  metrics_infos %>%
    filter(name %in% (models %>% unlist() %>% unique())) %T>%
    {pull(.,name) ->> metrics_names} %>% # Extract metrics names in the right order
    pull(path) %>%
    map(rast) %>%
    map(terra::project, epsg, method = "bilinear") %>%
    map(terra::resample, .[[1]]) %>%
    rast -> metrics

  names(metrics) <- metrics_names # Assign metrics names

  list(vect(paste0(wd, "/shp/roads.shp")),
       vect(paste0(wd, "/shp/waterbodies.shp"))) %>%
    map(project, epsg) -> masks

  # # Create subset area from drone lidar flight
  # st_read(paste0(wd, "/shp/placettes.shp")) %>%
  #   st_zm() %>%
  #   st_transform(epsg) %>%
  #   st_buffer(2000) %>%
  #   st_union() %>%
  #   vect() -> subset_area
  #
  # # Clip data for drone flight area
  # metrics %<>%
  #   crop(subset_area) %>%
  #   mask(subset_area)
  #
  # masks %<>%
  #   map(crop, subset_area) %>%
  #   map(mask, subset_area)

  # Perform segmentation
  for(m in models %>% names()){
    cat(paste0(m, "\n"))
    eFRI_segmentation(metrics = metrics[[models[[m]]]],
                      # masks = masks,
                      masks = NULL,
                      thresh = best_parameters_all %>% filter(forest == f, model_name == m) %>% pull(thresh),
                      spec = best_parameters_all %>% filter(forest == f, model_name == m) %>% pull(spec),
                      spat = best_parameters_all %>% filter(forest == f, model_name == m) %>% pull(spat),
                      method = "bs",
                      clean_nodata = TRUE,
                      output_path = paste0(wd, "/segmentations"),
                      output_name = m,
                      otb_dir = "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin") -> segmentation

    segmentation %>%
      st_write(paste0(wd, "/segmentations/data.gpkg"),
               layer = m,
               quiet = T)
  }
}
