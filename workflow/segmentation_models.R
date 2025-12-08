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





# # First ones (keep them?) ----
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
#
#
#
# Beatiful big segmentation story (BBSS) ----
# RMF
models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
               model_B = c("z_p95", "z_above2", "z_skew"),
               model_C = c("z_p95", "z_above2", "fractional_cover_2_5"),
               model_D = c("z_p95", "z_above2", "sagawi"),
               model_E = c("z_p95", "z_above2", "slope"),
               model_F = c("z_p95", "z_above2", "B6"),
               model_G = c("B6", "B8", "B11"),
               model_H = c("B2", "B3", "B4"),
               model_I = c("V_ha", "dens", "qmdbh"))

# OVF
models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
               model_B = c("z_p95", "z_above2", "z_skew"),
               model_D = c("z_p95", "z_above2", "sagawi"),
               model_E = c("z_p95", "z_above2", "slope"),
               model_F = c("z_p95", "z_above2", "B6"),
               model_G = c("B6", "B8", "B11"),
               model_H = c("B2", "B3", "B4"))





# Parameters ----
wd <- "D:/00_Ontario_eFRI/RMF" # wd
epsg <- "EPSG:2958" # RMF
#wd <- "D:/00_Ontario_eFRI/OVF" # wd
#epsg <- "EPSG:2959" # OVF
otb_dir <- "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin"
Sys.setenv(OTB_MAX_RAM_HINT = "65536") # 64 GO
Sys.setenv(OTB_MEMORY_AVAILABLE = "65536") # 64 GO
Sys.setenv(GDAL_CACHEMAX = "65536") # 64 GO





# Read metrics ----
list.files(paste0(wd, "/metrics"), full.names = TRUE) %>%
  read_metrics() -> metrics_infos

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

# Create subset area from drone lidar flight
st_read(paste0(wd, "/shp/placettes.shp")) %>%
  st_zm() %>%
  st_transform(epsg) %>%
  st_buffer(2000) %>%
  st_union() %>%
  vect() -> subset_area

# Clip data
metrics %<>%
  crop(subset_area) %>%
  mask(subset_area)

masks %<>%
  map(crop, subset_area) %>%
  map(mask, subset_area)





# Perform segmentation ----
map(names(models),
    function(x){
      cat(paste0(x, "\n"))
      eFRI_segmentation(metrics = metrics[[models[[x]]]],
                        # masks = masks,
                        masks = NULL,
                        #thresh = 47,
                        thresh = 70,
                        #spec = 0.6,
                        spec = 0.8,
                        #spat = 0.6,
                        spat = 0.4,
                        method = "bs",
                        clean_nodata = TRUE,
                        output_path = paste0(wd, "/segmentations/models_test_no_masks_cleaned_threshold"),
                        output_name = x,
                        otb_dir = "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin") -> segmentation

      segmentation %>%
        st_write(paste0(wd, "/segmentations/models_test_no_masks_cleaned_threshold/data.gpkg"),
                 layer = x,
                 quiet = T)
    })
