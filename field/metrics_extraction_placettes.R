# 🟡 Library ----
library(tidyverse)
library(terra)
library(sf)
library(lasR)
library(lidR)
library(lidRmetrics)
library(future)
library(eFRItools)





# 🟡 Setwd, parameters and functions ----
setwd("F:/Ontario_eFRI/02_donnees_traitees")
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/OLD_Rmarkdown/functions/fractional_transect.R")
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/OLD_Rmarkdown/functions/fractional_plot.R")
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/OLD_Rmarkdown/functions/fractional_cover.R")





# 🟡 Merge placettes ----
# RMF
map_dfr(1:5,
    function(x){
  st_read(paste0("./S", x ,"/placettes/01_corrige/placettes.shp"), quiet = T)
    }
) %>%
  st_write("placettes_RMF.shp")

# OVF
map_dfr(6:10,
        function(x){
          st_read(paste0("./S", x ,"/placettes/01_corrige/placettes.shp"), quiet = T)
        }
) %>%
  st_write("placettes_OVF.shp")





# 🟡 Computing ----
for(sub_sector in 1:10){
  cat(paste0("Sub-sector : S", sub_sector, "\n"))
  # Read data
  # Catalog
  readLAScatalog(paste0("./S", sub_sector ,"/nuages/02_norm_decimated")) -> ctg_normalized

  # Placettes
  st_read(paste0("./S", sub_sector ,"/placettes/01_corrige/placettes.shp"), quiet = T) -> placettes

  placette_data <- list()

  # Differents rayons
  for(radius_threshold in c(3.57, 5, 11.28, 14.10)){
    # Differents tri de zmin
    for(zmin_threshold in c(NA, 0.5, 1, 1.3, 2)){
      cat(paste0("Radius : ", radius_threshold, " / Zmin : ", zmin_threshold, "\n"))
      placettes %>%
        st_geometry() %>%
        plot_metrics(ctg_normalized,
                     ~lidRmetrics::metrics_set3(x = X, # Calculer toutes les metriques disponibles (package lidRmetrics)
                                                y = Y,
                                                z = Z,
                                                i = Intensity,
                                                ReturnNumber = ReturnNumber,
                                                NumberOfReturns = NumberOfReturns,
                                                zmin = zmin_threshold,
                                                threshold = c(2, 5),
                                                dz = 1,
                                                interval_count = 10,
                                                zintervals = c(0, 0.15, 2, 5, 10, 20, 30),
                                                pixel_size = 1,
                                                vox_size = 1,
                                                KeepReturns = c(1, 2, 3, 4)),
                     geometry = .,
                     radius = radius_threshold) %>%
        st_as_sf() %>%
        st_drop_geometry() %>%
        rename_with(.fn = ~ paste0(., "_radius_", radius_threshold*100, "_cm_zmin_", zmin_threshold*100, "_cm")) %>%
        append(placette_data) -> placette_data
    }

    # Tous les points sauf les points sols

    cat(paste0("Radius : ", radius_threshold, " / Zmin : drop ground \n"))

    opt_filter(ctg_normalized) <- "-drop_class 2"

    placettes %>%
      st_geometry() %>%
      plot_metrics(ctg_normalized,
                   ~lidRmetrics::metrics_set3(x = X, # Calculer toutes les metriques disponibles (package lidRmetrics)
                                              y = Y,
                                              z = Z,
                                              i = Intensity,
                                              ReturnNumber = ReturnNumber,
                                              NumberOfReturns = NumberOfReturns,
                                              zmin = NA,
                                              threshold = c(2, 5),
                                              dz = 1,
                                              interval_count = 10,
                                              zintervals = c(0, 0.15, 2, 5, 10, 20, 30),
                                              pixel_size = 1,
                                              vox_size = 1,
                                              KeepReturns = c(1, 2, 3, 4)),
                   geometry = .,
                   radius = radius_threshold) %>%
      st_as_sf() %>%
      st_drop_geometry() %>%
      rename_with(.fn = ~ paste0(., "_radius_", radius_threshold*100, "_cm_zmin_drop_ground_cm")) %>%
      append(placette_data) -> placette_data

    opt_filter(ctg_normalized) <- ""

    # Calculer le fractional cover sur un rayon de « radius_threshold »
    placettes %>%
      st_geometry() %>%
      st_as_sf() %>%
      fractional.plot(ctg_normalized,
                      radius = radius_threshold,
                      lower = c(1, 2, 3, 1),
                      upper = c(2, 3, 4, 4)) %>%
      st_drop_geometry() %>%
      rename_with(.fn = ~ paste0(., "_radius_", radius_threshold*100)) %>%
      append(placette_data) -> placette_data
  }
  # Calculer le fractional cover sur un transect Nord-Sud et Est-Ouest de 1m x 10m
  placettes %>%
    st_geometry() %>%
    st_as_sf() %>%
    fractional.transect(ctg_normalized,
                        length = 10,
                        width = 1,
                        lower = c(1, 2, 3, 1),
                        upper = c(2, 3, 4, 4)) %>%
    st_drop_geometry() %>%
    append(placette_data) -> placette_data

  placette_data %<>%
    bind_cols()

  placettes %>%
    bind_cols(placette_data) %>%
    saveRDS(paste0("./S", sub_sector ,"/placettes/01_corrige/placettes_data.rds"))
}





# Read test
map(1:10,
    ~{readRDS(paste0("./S", .x ,"/placettes/01_corrige/placettes_data.rds"))}) %>%
  map(dplyr::select) %>% # Maybe remove geometry if it doesn't work
  bind_rows() -> data

data %>%
  saveRDS("data.rds")

readRDS("data.rds") %>%
  dplyr::select(contains("drop_ground_cm")) %>%
  dplyr::select(contains("radius_1410_cm"))





# 🟡 Extract topographic metrics ----
map_dfr(c("RMF", "OVF"),
        function(sector){
          list.files(paste0("D:/00_Ontario_eFRI/", sector, "/metrics"), full.names = T) %>%
            read_metrics() %>%
            filter(name %in% c("aspect",
                               "dem",
                               "relative_topographic_position_100m",
                               "relative_topographic_position_200m",
                               "relative_topographic_position_400m",
                               "relative_topographic_position_60m",
                               "relative_topographic_position_800m",
                               "sagawi",
                               "slope",
                               "twi")) %T>%
            {pull(.,name) ->> metrics_names} %>% # Extract metrics names in the right order
            pull(path) %>%
            map(rast) %>%
            rast -> metrics

          names(metrics) <- metrics_names

          st_read(paste0("D:/00_Ontario_eFRI/", sector, "/shp/placettes_", sector, ".shp")) %>%
            st_zm() %>%
            st_buffer(11.28) %>%
            mutate_metrics(metrics, "mean") %>%
            st_drop_geometry()
        }) -> data_topo

data_topo %>%
  saveRDS("D:/00_Ontario_eFRI/data_topo.rds")
