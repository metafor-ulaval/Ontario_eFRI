library(tidyverse)
library(magrittr)
library(sf)
library(rmapshaper)
library(terra)
# pak::pak("metafor-ulaval/eFRItools")
library(eFRItools)



# # Merge zq95, zq90 into one raster per forest for final statistics ----
# dfa <- st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp")
#
# for(f in c("RMF", "OVF")){
#
#   cat(f, "\n")
#
#   dfa %>%
#     filter(forest == f) %>%
#     pull(sector) %>%
#     map(function(x){paste0("D:/00_Ontario_eFRI/random_forest/metriques/", x, "/metrics/zq95.tif")}) %>%
#     map(rast) %>%
#     sprc() %>%
#     terra::merge() -> zq95
#
#   names(zq95) <- "zq95"
#
#   zq95 %>%
#     writeRaster(paste0("./results/", f, "_zq95_prediction.tif"))
#
# }
#
# for(f in c("RMF", "OVF")){
#
#   cat(f, "\n")
#
#   dfa %>%
#     filter(forest == f) %>%
#     pull(sector) %>%
#     map(function(x){paste0("D:/00_Ontario_eFRI/random_forest/metriques/", x, "/metrics/zq90.tif")}) %>%
#     map(rast) %>%
#     sprc() %>%
#     terra::merge() -> zq90
#
#   names(zq90) <- "zq90"
#
#   zq90 %>%
#     writeRaster(paste0("./results/", f, "_zq90_prediction.tif"))
#
# }
#
#
#
#

# Final statistics ----
st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp") %>%
  filter(sector != "S5") -> dfa

# x <- "RMF"
# y <- "zq90"
# z <- "model_B"
# s <- 5000

map_dfr(dfa %>% pull(forest) %>% unique,
        function(x){

          st_read(paste0("D:/00_Ontario_eFRI/", x, "/shp/PolygonForest.shp"), quiet = T) %>%
            st_cast("MULTIPOLYGON") %>%
            st_cast("POLYGON") -> pf

          st_layers(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg")) %>%
            as_tibble %>%
            pull(name) -> models

          dfa %>%
            filter(forest == x) -> dfa_temp

          map_dfr(c("vmb_ha", "st", "dens", "dhpq", "zq90"),
                  function(y){

                    rast(paste0("D:/00_Ontario_eFRI/random_forest/results/", x, "_", y, "_prediction.tif")) -> metric

                    st_crs(metric) -> ref_crs

                    dfa_temp %<>%
                      st_transform(ref_crs)

                    metric %<>%
                      crop(dfa_temp)

                    ifel(!is.na(metric), 1, NA) %>%
                      as.polygons() %>%
                      st_as_sf() -> metric_mask

                    map_dfr(models,
                            function(z){

                              cat(x, y, z, "\n")

                              st_read(paste0("D:/00_Ontario_eFRI/", x, "/segmentations/data.gpkg"), layer = z, quiet = T) %>%
                                st_cast("MULTIPOLYGON") %>%
                                st_cast("POLYGON") -> segmentation

                              pf %>%
                                st_filter(metric_mask %>% st_transform(st_crs(pf))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(metric_mask) %>%
                                mutate(area_clip = as.numeric(st_area(.))) %>%
                                dplyr::select(area, area_clip) %>%
                                mutate(prop = as.numeric(area_clip/area*100)) %>%
                                mutate(sd = exactextractr::exact_extract(metric, ., fun = "stdev")) -> pf_temp

                              segmentation %>%
                                st_filter(metric_mask %>% st_transform(st_crs(segmentation))) %>%
                                st_transform(ref_crs) %>%
                                mutate(area = as.numeric(st_area(.))) %>%
                                ms_clip(metric_mask) %>%
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

                                        proportion_variance_explained(pf_temp_sup, metric) -> pve_forest
                                        proportion_variance_explained(segmentation_temp_sup, metric) -> pve_segmentation

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
  write.csv("D:/00_Ontario_eFRI/segmentation/results/results_final.csv")

results %>%
  saveRDS("D:/00_Ontario_eFRI/segmentation/results/results_final.rds")
