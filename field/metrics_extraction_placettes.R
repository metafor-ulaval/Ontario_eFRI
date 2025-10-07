# 🟡 Library ----
library(tidyverse)
library(terra)
library(sf)
library(lasR)
library(lidR)
library(lidRmetrics)
library(future)


# 🟡 Setwd, parameters and functions ----
setwd("F:/Ontario_eFRI/02_donnees_traitees")
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/OLD_Rmarkdown/functions/fractional_transect.R")
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/OLD_Rmarkdown/functions/fractional_plot.R")
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/OLD_Rmarkdown/functions/fractional_cover.R")

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
        st_drop_geometry() %>%
        dplyr::select(-Name) %>%
        rename_with(.fn = ~ paste0(., "_radius_", radius_threshold*100, "_cm_zmin_", zmin_threshold*100, "_cm")) %>%
        append(placette_data) -> placette_data
    }

    # Tous les points sauf les points sols

    cat(paste0("Radius : ", radius_threshold, " / Zmin : drop ground \n"))

    opt_filter(ctg_normalized) <- "-drop_class 2"

    placettes %>%
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
      st_drop_geometry() %>%
      dplyr::select(-Name) %>%
      rename_with(.fn = ~ paste0(., "_radius_", radius_threshold*100, "_cm_zmin_drop_ground_cm")) %>%
      append(placette_data) -> placette_data

    opt_filter(ctg_normalized) <- ""

    # Calculer le fractional cover sur un rayon de « radius_threshold »
    placettes %>%
      fractional.plot(ctg_normalized,
                      radius = radius_threshold,
                      lower = c(1, 2, 3, 1),
                      upper = c(2, 3, 4, 4)) %>%
      st_drop_geometry() %>%
      dplyr::select(-Name) %>%
      rename_with(.fn = ~ paste0(., "_radius_", radius_threshold*100)) %>%
      append(placette_data) -> placette_data
  }
  # Calculer le fractional cover sur un transect Nord-Sud et Est-Ouest de 1m x 10m
  placettes %>%
    fractional.transect(ctg_normalized,
                        length = 10,
                        width = 1,
                        lower = c(1, 2, 3, 1),
                        upper = c(2, 3, 4, 4)) %>%
    append(placette_data) -> placette_data

  placette_data %>%
    bind_cols() %>%
    st_write(paste0("./", sub_sector ,"/placettes/01_corrige/placettes_data.gpkg"), "placettes_data", quiet = T)
}
