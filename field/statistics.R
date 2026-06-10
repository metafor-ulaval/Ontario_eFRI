library(tidyverse)
library(magrittr)
library(sf)
library(rmapshaper)
library(terra)
# pak::pak("metafor-ulaval/eFRItools")
library(eFRItools)




# Internal homogeneity sensibility analysis ----
dfa <- st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp")

x <- "RMF"
y <- "vmbha"
z <- "model_B"
s <- 2000

map_dfr(dfa %>% pull(forest) %>% unique,
        function(x){

          pf <- st_read(paste0("D:/00_Ontario_eFRI/", x, "/shp/PolygonForest.shp"), quiet = T) %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")
          models <- st_layers(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg")) %>% as_tibble %>% pull(name)
          dfa_temp <- dfa %>% filter(forest == x)

          map_dfr(c("vmbha", "st", "dens", "dhpq"),
                  function(y){

                    metric <- rast(paste0("D:/dossier_remise_fin_contrat_PY/ontario/predictions/2026_05_14/", y, "_", x, ".tif"))

                    map_dfr(models,
                            function(z){

                              cat(x, y, z, "\n")

                              segmentation <- st_read(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg"), layer = z, quiet = T) %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")

                              ref_crs <- st_crs(segmentation)
                              dfa_temp %<>%
                                st_transform(ref_crs)

                              pf %>%
                                st_filter(dfa_temp %>% st_transform(st_crs(pf))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(dfa_temp) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) %>%
                                mutate(sd = exactextractr::exact_extract(metric, ., fun = "stdev")) %>%
                                drop_na(sd) -> pf_temp

                              segmentation %>%
                                st_filter(dfa_temp %>% st_transform(st_crs(segmentation))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(dfa_temp) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) %>%
                                mutate(sd = exactextractr::exact_extract(metric, ., fun = "stdev")) %>%
                                drop_na(sd) -> segmentation_temp

                              # ggplot() +
                              #   geom_sf(data = pf_temp, col = "green", fill = NA) +
                              #   geom_sf(data = segmentation_temp, col = "blue", fill = NA) +
                              #   geom_sf(data = dfa_temp, col = "red", fill = NA)

                              map_dfr(1000*2^seq(1, 5, 1),
                                      function(s){

                                        cat(s, "- ")

                                        pf_temp %>%
                                          filter(area_clip >= s) -> pf_temp_sup

                                        segmentation_temp %>%
                                          filter(area_clip >= s) -> segmentation_temp_sup

                                        tibble(forest = x,
                                               metric = y,
                                               model = z,
                                               sup_min = s,
                                               n_forest = nrow(pf_temp_sup),
                                               n_segmentation = nrow(segmentation_temp_sup),
                                               ih_forest = sum(pf_temp_sup$area_clip*pf_temp_sup$sd)/sum(pf_temp_sup$area_clip),
                                               ih_segmentation = sum(segmentation_temp_sup$area_clip*segmentation_temp_sup$sd)/sum(segmentation_temp_sup$area_clip)) -> results

                                        return(results)

                                      }) -> results

                              cat("\n")

                              return(results)

                            }) -> results

                    return(results)

                  }) -> results

          return(results)

        }) -> results

results %>%
  saveRDS("D:/00_Ontario_eFRI/segmentation/results/results_internal_homogeneity.rds")





# Proportion of variance explained sensibility analysis ----
dfa <- st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp")

x <- "RMF"
y <- "vmbha"
z <- "model_A"
s <- 32000

map_dfr(dfa %>% pull(forest) %>% unique,
        function(x){

          pf <- st_read(paste0("D:/00_Ontario_eFRI/", x, "/shp/PolygonForest.shp"), quiet = T) %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")
          models <- st_layers(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg")) %>% as_tibble %>% pull(name)
          dfa_temp <- dfa %>% filter(forest == x)

          map_dfr(c("vmbha", "st", "dens", "dhpq"),
                  function(y){

                    metric <- rast(paste0("D:/dossier_remise_fin_contrat_PY/ontario/predictions/2026_05_14/", y, "_", x, ".tif"))

                    map_dfr(models,
                            function(z){

                              cat(x, y, z, "\n")

                              segmentation <- st_read(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg"), layer = z, quiet = T) %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")

                              ref_crs <- st_crs(segmentation)
                              dfa_temp %<>%
                                st_transform(ref_crs)

                              pf %>%
                                st_filter(dfa_temp %>% st_transform(st_crs(pf))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(dfa_temp) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) -> pf_temp

                              segmentation %>%
                                st_filter(dfa_temp %>% st_transform(st_crs(segmentation))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(dfa_temp) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) -> segmentation_temp

                              # ggplot() +
                              #   geom_sf(data = pf_temp, col = "green", fill = NA) +
                              #   geom_sf(data = segmentation_temp, col = "blue", fill = NA) +
                              #   geom_sf(data = dfa_temp, col = "red", fill = NA)

                              map_dfr(1000*2^seq(1, 5, 1),
                                      function(s){

                                        cat(s, "- ")

                                        pf_temp %>%
                                          filter(area_clip >= s) -> pf_temp_sup

                                        segmentation_temp %>%
                                          filter(area_clip >= s) -> segmentation_temp_sup

                                        pve_forest <- proportion_variance_explained(pf_temp_sup, metric)
                                        pve_segmentation <- proportion_variance_explained(segmentation_temp_sup, metric)

                                        tibble(forest = x,
                                               metric = y,
                                               model = z,
                                               sup_min = s,
                                               n_forest = nrow(pf_temp_sup),
                                               n_segmentation = nrow(segmentation_temp_sup),
                                               pve_forest = pve_forest,
                                               pve_segmentation = pve_segmentation) -> results

                                        return(results)

                                      }) -> results

                              cat("\n")

                              return(results)

                            }) -> results

                    return(results)

                  }) -> results

          return(results)

        }) -> results

results %>%
  saveRDS("D:/00_Ontario_eFRI/segmentation/results/results_proportion_variance_explained.rds")





# Final statistics ----
dfa <- st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp")

# x <- "RMF"
# y <- "vmbha"
# z <- "model_B"
# s <- 5000

map_dfr(dfa %>% pull(forest) %>% unique,
        function(x){

          pf <- st_read(paste0("D:/00_Ontario_eFRI/", x, "/shp/PolygonForest.shp"), quiet = T) %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")
          models <- st_layers(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg")) %>% as_tibble %>% pull(name)
          dfa_temp <- dfa %>% filter(forest == x)

          map_dfr(c("vmbha", "st", "dens", "dhpq"),
                  function(y){

                    metric <- rast(paste0("D:/dossier_remise_fin_contrat_PY/ontario/predictions/2026_05_14/", y, "_", x, ".tif"))

                    map_dfr(models,
                            function(z){

                              cat(x, y, z, "\n")

                              segmentation <- st_read(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg"), layer = z, quiet = T) %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")

                              ref_crs <- st_crs(segmentation)
                              dfa_temp %<>%
                                st_transform(ref_crs)

                              pf %>%
                                st_filter(dfa_temp %>% st_transform(st_crs(pf))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(dfa_temp) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) %>%
                                mutate(sd = exactextractr::exact_extract(metric, ., fun = "stdev")) -> pf_temp

                              segmentation %>%
                                st_filter(dfa_temp %>% st_transform(st_crs(segmentation))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(dfa_temp) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) %>%
                                mutate(sd = exactextractr::exact_extract(metric, ., fun = "stdev")) -> segmentation_temp

                              # ggplot() +
                              #   geom_sf(data = pf_temp, col = "green", fill = NA) +
                              #   geom_sf(data = segmentation_temp, col = "blue", fill = NA) +
                              #   geom_sf(data = dfa_temp, col = "red", fill = NA)

                              map_dfr(c(5000, 10000, 15000),
                                      function(s){

                                        cat(s, "- ")

                                        pf_temp %>%
                                          filter(area_clip >= s) -> pf_temp_sup

                                        segmentation_temp %>%
                                          filter(area_clip >= s) -> segmentation_temp_sup

                                        pve_forest <- proportion_variance_explained(pf_temp_sup, metric)
                                        pve_segmentation <- proportion_variance_explained(segmentation_temp_sup, metric)

                                        tibble(forest = x,
                                               metric = y,
                                               model = z,
                                               sup_min = s,
                                               n_forest = nrow(pf_temp_sup),
                                               n_segmentation = nrow(segmentation_temp_sup),
                                               pve_forest = pve_forest,
                                               pve_segmentation = pve_segmentation,
                                               ih_forest = sum(pf_temp_sup$area_clip*pf_temp_sup$sd)/sum(pf_temp_sup$area_clip),
                                               ih_segmentation = sum(segmentation_temp_sup$area_clip*segmentation_temp_sup$sd)/sum(segmentation_temp_sup$area_clip)) -> results

                                        return(results)

                                      }) -> results

                              cat("\n")

                              return(results)

                            }) -> results

                    return(results)

                  }) -> results

          return(results)

        }) -> results

results %>%
  saveRDS("D:/00_Ontario_eFRI/segmentation/results/results_final.rds")
