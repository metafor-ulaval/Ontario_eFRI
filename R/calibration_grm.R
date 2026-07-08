# 🟡 Library ----
library(tidyverse)
library(future.apply)
library(future.callr)
library(future)
library(magrittr)
library(terra)
library(rmapshaper)
library(sf)
library(lasR)
library(lidR)
library(lidRmetrics)
library(plotly)
library(htmlwidgets)
library(eFRItools)
library(basemapR)
library(patchwork)

sample_farthest_point <- function(nb,
                                  res,
                                  sampling_area,
                                  seed) {

  # ------------------------------------------------------------
  # 1. Validate inputs
  # ------------------------------------------------------------
  if (!inherits(sampling_area, "sf")) {
    stop("Sampling_area must be sf objects.")
  }

  # ------------------------------------------------------------
  # 2. Prepare a permanent empty raster template (created ONCE)
  # ------------------------------------------------------------
  bbox <- sf::st_bbox(sampling_area)

  rast_template <- terra::rast(
    xmin  = bbox["xmin"],
    xmax  = bbox["xmax"],
    ymin  = bbox["ymin"],
    ymax  = bbox["ymax"],
    resolution = res,
    crs = sf::st_crs(sampling_area)$wkt
  )

  terra::values(rast_template) <- 1

  # ------------------------------------------------------------
  # 3. Select the first point randomly
  # ------------------------------------------------------------
  set.seed(seed)

  pt <- st_sample(sampling_area, 1, type = "random")

  # ------------------------------------------------------------
  # 4. Iteratively select farthest points
  # ------------------------------------------------------------
  for (i in 1:(nb-1)) {

    # Compute distance raster FROM currently selected points
    dist_raster <- terra::distance(rast_template, terra::vect(pt))

    # Mask raster
    dist_raster_masked <- terra::mask(dist_raster, sampling_area)

    # Plot results
    plot(dist_raster_masked)
    plot(pt, add = T, col = "red", pch = 16)

    # Locate farthest cell location within raster
    max_cell <- which.max(values(dist_raster_masked))
    farthest_cell <- rast(dist_raster_masked)
    values(farthest_cell) <- NA
    values(farthest_cell)[max_cell] <- 1

    # Randomly sample one point inside this cell
    farthest_cell_poly <- terra::as.polygons(farthest_cell)
    farthest_cell_poly <- sf::st_as_sf(farthest_cell_poly)
    set.seed(seed)
    new_pt <- st_sample(farthest_cell_poly, 1, type = "random")

    # Plot results
    plot(new_pt, add = T, col = "blue", pch = 16)

    pt <- c(pt, new_pt)

  }

  return(pt)
}






# 🟡 Generate base data for segmentation ----
f <- "RMF"
m <- c("z_p95", "z_above2", "z_cv", "z_skew", "sagawi", "B6", "B8", "B11")
inner_radius <- sqrt(10000^2*2)/2

for(f in c("RMF", "OVF")){

  st_read(paste0("D:/00_Ontario_eFRI/", f, "/shp/PolygonForest.shp"), quiet = T) -> polygon_forest

  st_read(paste0("D:/00_Ontario_eFRI/", f, "/ctg/ctg.shp"), quiet = T) %>%
    st_transform(st_crs(polygon_forest)) -> ctg

  ctg %>%
    st_geometry() %>%
    st_buffer(1000) %>%
    st_union() %>%
    st_buffer(-1000) -> area

  polygon_forest %>%
    filter(POLYTYPE == "FOR") %>%
    st_geometry() %>%
    st_as_sf() %>%
    sf::st_cast("MULTIPOLYGON") %>%
    st_cast("POLYGON") %>%
    st_filter(area) %>%
    st_buffer(500) %>%
    st_union() %>%
    st_buffer(-500) %>%
    st_buffer(-inner_radius)  %>%
    st_as_sf() -> sampling_inner_area

  sample_farthest_point(nb = 10,
                        res = 100,
                        sampling_area = sampling_inner_area,
                        seed = 1) %>%
    st_as_sf() -> sampling_point

  map_dfr(1:10,
          function(x){
            sampling_point[x,] %>%
              st_buffer(2500) %>%
              st_as_sf() %>%
              st_bbox() %>%
              st_as_sfc() %>%
              st_as_sf()
          }) %>%
    st_as_sf() -> sampling_square

  sampling_square %>%
    st_transform(st_crs(polygon_forest)) %>%
    st_write(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/sampling_square_", f, ".shp"))

  polygon_forest %>%
    st_intersection(sampling_square) %>%
    rowid_to_column("id") -> polygon_forest

  polygon_forest %>%
    st_write(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/polygon_forest_", f,".shp"))

  polygon_forest %>%
    sf::st_cast("MULTIPOLYGON") %>%
    st_cast("POLYGON") %>%
    filter(POLYTYPE == "FOR") -> polygon_forest_FOR

  polygon_forest %>%
    filter(!id %in% polygon_forest_FOR$id) %>%
    st_geometry() %>%
    sf::st_cast("MULTIPOLYGON") %>%
    st_cast("POLYGON") %>%
    st_as_sf() -> polygon_forest_NON_FOR

  list.files(paste0("D:/00_Ontario_eFRI/", f, "/metrics"), full.names = TRUE) %>%
    read_metrics() %>%
    filter(name %in% m) %T>%
    {pull(.,name) ->> metrics_names} %>% # Extract metrics names in the right order
    pull(path) %>%
    map(rast) %>%
    map(terra::project, paste0("EPSG:", st_crs(polygon_forest)$epsg), method = "bilinear") %>%
    map(terra::resample, .[[1]]) %>%
    rast -> metrics

  names(metrics) <- metrics_names

  metrics %<>%
    terra::crop(polygon_forest_FOR) %>%
    terra::mask(polygon_forest_FOR)

  # Save objects
  metrics %>%
    writeRaster(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/metrics_", f,".tif"))

  polygon_forest_FOR %>%
    st_write(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/polygon_forest_FOR_", f,".shp"))

  polygon_forest_NON_FOR %>%
    st_write(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/masks_", f,".shp"))

}





# 🟡 Segmentation in parallel ----
plan(multisession, workers = round(ncores()*0.25))

forest <- c("RMF", "OVF")
thresh <- seq(10, 100, 10)
spec <- seq(0.1, 1, 0.1)
spat <- seq(0.1, 1, 0.1)
models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
               model_B = c("z_p95", "z_above2", "z_skew"),
               model_D = c("z_p95", "z_above2", "sagawi"),
               model_F = c("z_p95", "z_above2", "B6"),
               model_G = c("B6", "B8", "B11"))

models_df <- data.frame(model_name = names(models),
                        model_var  = sapply(models, paste, collapse = ", "),
                        stringsAsFactors = FALSE)

params <- merge(expand.grid(forest = forest,
                            thresh = thresh,
                            spec = spec,
                            spat = spat,
                            stringsAsFactors = FALSE),
                models_df)

# Remove done files
params %>%
  mutate(done = file.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", forest, "/", model_name, "/segmentation_thresh", thresh, "_spec", spec*100, "_spat", spat*100, "_smoothed.shp"))) %>%
  filter(!done) -> params

# Remove remaining raster files
list.files(path = "D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation", full.names = T) %>%
  map(list.files, full.names = T) %>%
  unlist() %>%
  map(list.files, pattern = "\\.tif$", full.names = T) %>%
  unlist() %>%
  file.remove()

# params %>%
#   filter(forest == "RMF" & model_name == "model_F" & thresh == 60 & spec == 0.2 & spat == 0.1) -> param

# p <- 889
# x <- 1

future_lapply(seq_len(nrow(params)), function(p){

  library(eFRItools)
  library(sf)
  library(terra)
  library(stringr)

  param = params[p,]

  cat("Forest :", param$forest, "/ Model :", param$model_name, "/ Threshold :", param$thresh, "/ Spec :", param$spec, "/ Spat :", param$spat, "\n")

  if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest))){
    dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest))
  }

  if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name))){
    dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name))
  }

  output_name <- paste0("_thresh", param$thresh,
                        "_spec", param$spec*100,
                        "_spat", param$spat*100)

  st_read(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/sampling_square_", param$forest, ".shp"),
          quiet = T) -> sampling_square

  rast(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/metrics_", param$forest,".tif")) -> metrics

  st_read(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/masks_", param$forest,".shp"),
          quiet = T) -> masks

  map_dfr(seq_len(nrow(sampling_square)),
          function(x){

            metrics[[strsplit(param$model_var, ", ")[[1]]]] %>%
              terra::crop(sampling_square[x,]) %>%
              terra::mask(sampling_square[x,]) -> metrics_temp

            masks %>%
              sf::st_filter(sampling_square[x,]) -> masks_temp

            eFRI_segmentation(metrics = metrics_temp,
                              masks = list(vect(masks_temp)),
                              thresh = param$thresh,
                              spec = param$spec,
                              spat = param$spat,
                              method = "bs",
                              clean_nodata = TRUE,
                              output_path = paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name),
                              output_name = paste0("segmentation", output_name, "_it", x),
                              otb_dir = "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin") %>%
              sf::st_geometry() %>%
              sf::st_as_sf() -> segmentation

            sf::st_geometry(segmentation) <- "geometry"

            return(segmentation)

          }) -> segmentation_merged

  list.files(path = paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name),
             pattern = "\\.tif$",
             full.names = T) %>%
    stringr::str_subset(paste0("segmentation", output_name)) %>%
    file.remove()

  st_write(segmentation_merged,
           paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name, "/segmentation", output_name, ".shp"),
           append = FALSE,
           quiet = T)

  chunk_size <- 1000
  segmentation_chunks <- split(segmentation_merged,
                               ceiling(seq_len(nrow(segmentation_merged)) / chunk_size))

  map_dfr(seq_along(segmentation_chunks),
          function(chunk_no) {

            cat("Processing chunk", chunk_no, "out of", length(segmentation_chunks), "\n")

            chunk <- segmentation_chunks[[chunk_no]]

            chunk %>%
              smoothr::smooth(method = "ksmooth", smoothness = 1) %>%
              sf::st_make_valid() -> chunk_smoothed

              return(chunk_smoothed)

          }) -> segmentation_smoothed

  segmentation_smoothed %>%
    st_cast("MULTIPOLYGON", warn = FALSE) %>%
    st_make_valid() -> segmentation_smoothed

  segmentation_smoothed[st_geometry_type(segmentation_smoothed) %in% c("POLYGON", "MULTIPOLYGON"),] -> segmentation_smoothed

  st_write(segmentation_smoothed,
           paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name, "/segmentation", output_name, "_smoothed.shp"),
           append = FALSE,
           quiet = T)

  return(TRUE)

})





# 🟡 Compute results in parallel ----
forest <- c("RMF", "OVF")
thresh <- seq(10, 100, 10)
spec <- seq(0.1, 1, 0.1)
spat <- seq(0.1, 1, 0.1)
models <- list(model_A = c("z_p95", "z_above2", "z_cv"),
               model_B = c("z_p95", "z_above2", "z_skew"),
               model_D = c("z_p95", "z_above2", "sagawi"),
               model_F = c("z_p95", "z_above2", "B6"),
               model_G = c("B6", "B8", "B11"))

models_df <- data.frame(model_name = names(models),
                        model_var  = sapply(models, paste, collapse = ", "),
                        stringsAsFactors = FALSE)

params <- merge(expand.grid(forest = forest,
                            thresh = thresh,
                            spec = spec,
                            spat = spat,
                            stringsAsFactors = FALSE),
                models_df)

# Remove done files
params %>%
  mutate(done = file.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", forest, "/", model_name, "/tables/results_thresh", thresh, "_spec", spec*100, "_spat", spat*100, ".rds"))) %>%
  filter(!done) -> params

# p <- 1

params %>%
  filter(forest == "OVF" & model_name == "model_F" & thresh == 70 & spec == 0.6 & spat == 0.9) -> param

batch_size <- ncores()  # petit !

indices <- split(seq_len(nrow(params)),
                 ceiling(seq_along(seq_len(nrow(params))) / batch_size))

for(batch in indices){

  plan(callr, workers = round(lasR::ncores()*0.25))

  future.apply::future_lapply(batch, function(p){

    tryCatch({

      library(tidyverse)
      library(sf)
      library(terra)
      library(overlapping)
      library(patchwork)
      library(rmapshaper)
      library(magrittr)

      sspe <- function(actual, predicted) {
        (predicted - actual) / ((abs(actual) + abs(predicted)) / 2)
      }

      param = params[p,]

      cat("Forest :", param$forest, "/ Model :", param$model_name, "/ Threshold :", param$thresh, "/ Spec :", param$spec, "/ Spat :", param$spat, "\n")

      if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest))){
        dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest))
      }

      if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name))){
        dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name))
      }

      if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_nam, "/tables"))){
        dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/tables"))
      }

      if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/figures"))){
        dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/figures"))
      }

      if(!dir.exists(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/maps"))){
        dir.create(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/maps"))
      }

      output_name <- paste0("_thresh", param$thresh,
                            "_spec", param$spec*100,
                            "_spat", param$spat*100)

      terra::rast(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/metrics_", param$forest, ".tif")) -> metrics

      sf::st_read(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/base_data/polygon_forest_FOR_", param$forest, ".shp"), quiet = T) %>%
        dplyr::select() %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON") %>%
        mutate(area = as.numeric(sf::st_area(.)),
               perimeter = as.numeric(sf::st_perimeter(.)),
               perimeter_area_ratio = perimeter/area,
               shape_index = perimeter / sqrt(pi * area)) %>%
        filter_if(is.numeric, is.finite) -> polygon_forest

      sf::st_read(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name, "/segmentation", output_name, ".shp"), quiet = T) %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON") %>%
        mutate(area = as.numeric(st_area(.)),
               perimeter = as.numeric(st_perimeter(.)),
               perimeter_area_ratio = perimeter/area,
               shape_index = perimeter / sqrt(pi * area)) %>%
        filter_if(is.numeric, is.finite) -> segmentation

      sf::st_read(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/segmentation/", param$forest, "/", param$model_name, "/segmentation", output_name, "_smoothed.shp"), quiet = T) %>%
        sf::st_make_valid() %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON") %>%
        mutate(area = as.numeric(st_area(.)),
               perimeter = as.numeric(st_perimeter(.)),
               perimeter_area_ratio = perimeter/area,
               shape_index = perimeter / sqrt(pi * area)) %>%
        filter_if(is.numeric, is.finite) -> segmentation_smoothed

      polygon_forest %>%
        dplyr::select(area) %>%
        terra::rasterize(metrics, field = "area") %>%
        terra::values(na.rm = TRUE) -> polygon_forest_area_raster

      polygon_forest_area_raster[is.finite(polygon_forest_area_raster)] -> polygon_forest_area_raster

      segmentation %>%
        dplyr::select(area) %>%
        terra::rasterize(metrics, field = "area") %>%
        terra::values(na.rm = TRUE) -> segmentation_area_raster

      segmentation_area_raster[is.finite(segmentation_area_raster)] -> segmentation_area_raster

      overlap_area_raster_OV <- overlapping::overlap(list(segmentation_area_raster, polygon_forest_area_raster))$OV

      # ggplot() +
      #   geom_density(data = segmentation, aes(x = area/10000), color = "#D81B60") +
      #   geom_density(data = polygon_forest, aes(x = area/10000), color = "#0B4D00") +
      #   theme_bw() -> p_area_poly

      ggplot() +
        geom_histogram(data = segmentation %>% dplyr::filter(is.finite(area)),
                       aes(x = area/10000,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#D81B60",
                       alpha = 0.3) +
        geom_density(data = segmentation %>% dplyr::filter(is.finite(area)),
                     aes(x = area/10000),
                     color = "#D81B60",
                     linewidth = 1) +

        geom_histogram(data = polygon_forest %>% dplyr::filter(is.finite(area)),
                       aes(x = area/10000,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#0B4D00",
                       alpha = 0.3) +
        geom_density(data = polygon_forest %>% dplyr::filter(is.finite(area)),
                     aes(x = area/10000),
                     color = "#0B4D00",
                     linewidth = 1) +
        theme_bw() -> p_area_poly

      # ggplot() +
      #   geom_density(data = segmentation_area_raster, aes(x = area/10000), color = "#D81B60") +
      #   geom_density(data = polygon_forest_area_raster, aes(x = area/10000), color = "#0B4D00") +
      #   theme_bw() -> p_area_raster

      ggplot() +
        geom_histogram(aes(x = segmentation_area_raster/10000,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#D81B60",
                       alpha = 0.3) +
        geom_density(aes(x = segmentation_area_raster/10000),
                     color = "#D81B60",
                     linewidth = 1) +

        geom_histogram(aes(x = polygon_forest_area_raster/10000,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#0B4D00",
                       alpha = 0.3) +
        geom_density(aes(x = polygon_forest_area_raster/10000),
                     color = "#0B4D00",
                     linewidth = 1) +
        theme_bw() -> p_area_raster

      p_area_poly / p_area_raster -> p_area

      ggsave(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/figures/area", output_name, ".jpg"),
             plot = p_area,
             width = 15,
             height = 10,
             dpi = 150)

      graphics.off()
      rm(polygon_forest_area_raster,
         segmentation_area_raster,
         p_area_poly,
         p_area_raster,
         p_area)
      gc()

      polygon_forest %>%
        dplyr::select(shape_index) %>%
        terra::rasterize(metrics, field = "shape_index") %>%
        terra::values(na.rm = TRUE) -> polygon_forest_shape_index_raster

      polygon_forest_shape_index_raster[is.finite(polygon_forest_shape_index_raster)] -> polygon_forest_shape_index_raster

      segmentation %>%
        dplyr::select(shape_index) %>%
        terra::rasterize(metrics, field = "shape_index") %>%
        terra::values(na.rm = TRUE) -> segmentation_shape_index_raster

      segmentation_shape_index_raster[is.finite(segmentation_shape_index_raster)] -> segmentation_shape_index_raster

      segmentation_smoothed %>%
        dplyr::select(shape_index) %>%
        terra::rasterize(metrics, field = "shape_index") %>%
        terra::values(na.rm = TRUE) -> segmentation_smoothed_shape_index_raster

      segmentation_smoothed_shape_index_raster[is.finite(segmentation_smoothed_shape_index_raster)] -> segmentation_smoothed_shape_index_raster

      overlap_shape_raster_OV <- overlapping::overlap(list(segmentation_shape_index_raster, polygon_forest_shape_index_raster))$OV

      overlap_shape_smoothed_raster_OV <- overlapping::overlap(list(segmentation_smoothed_shape_index_raster, polygon_forest_shape_index_raster))$OV

      # ggplot() +
      #   geom_density(data = segmentation, aes(x = shape_index), color = "#D81B60") +
      #   geom_density(data = segmentation_smoothed, aes(x = shape_index), color = "#1E88E5") +
      #   geom_density(data = polygon_forest, aes(x = shape_index), color = "#0B4D00") +
      #   theme_void() -> p_shape_poly

      ggplot() +
        geom_histogram(data = segmentation %>% dplyr::filter(is.finite(shape_index)),
                       aes(x = shape_index,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#D81B60",
                       alpha = 0.3) +
        geom_density(data = segmentation %>% dplyr::filter(is.finite(shape_index)),
                     aes(x = shape_index),
                     color = "#D81B60",
                     linewidth = 1) +

        geom_histogram(data = segmentation_smoothed %>% dplyr::filter(is.finite(shape_index)),
                       aes(x = shape_index,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#1E88E5",
                       alpha = 0.3) +
        geom_density(data = segmentation_smoothed %>% dplyr::filter(is.finite(shape_index)),
                     aes(x = shape_index),
                     color = "#1E88E5",
                     linewidth = 1) +

        geom_histogram(data = polygon_forest %>% dplyr::filter(is.finite(shape_index)),
                       aes(x = shape_index,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#0B4D00",
                       alpha = 0.3) +
        geom_density(data = polygon_forest %>% dplyr::filter(is.finite(shape_index)),
                     aes(x = shape_index),
                     color = "#0B4D00",
                     linewidth = 1) +
        theme_bw() -> p_shape_poly

      # ggplot() +
      #   geom_density(data = segmentation_shape_index_raster, aes(x = shape_index), color = "#D81B60") +
      #   geom_density(data = segmentation_smoothed_shape_index_raster, aes(x = shape_index), color = "#1E88E5") +
      #   geom_density(data = polygon_forest_shape_index_raster, aes(x = shape_index), color = "#0B4D00") +
      #   theme_void() -> p_shape_raster

      ggplot() +
        geom_histogram(aes(x = segmentation_shape_index_raster,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#D81B60",
                       alpha = 0.3) +
        geom_density(aes(x = segmentation_shape_index_raster),
                     color = "#D81B60",
                     linewidth = 1) +

        geom_histogram(aes(x = segmentation_smoothed_shape_index_raster,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#1E88E5",
                       alpha = 0.3) +
        geom_density(aes(x = segmentation_smoothed_shape_index_raster),
                     color = "#1E88E5",
                     linewidth = 1) +

        geom_histogram(aes(x = polygon_forest_shape_index_raster,
                           y = after_stat(density)),
                       bins = 100,
                       fill = "#0B4D00",
                       alpha = 0.3) +
        geom_density(aes(x = polygon_forest_shape_index_raster),
                     color = "#0B4D00",
                     linewidth = 1) +
        theme_bw() -> p_shape_raster

      p_shape_poly / p_shape_raster -> p_shape

      ggsave(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/figures/shape", output_name, ".jpg"),
             plot = p_shape,
             width = 15,
             height = 10,
             dpi = 150)

      graphics.off()
      rm(polygon_forest_shape_index_raster,
         segmentation_shape_index_raster,
         segmentation_smoothed_shape_index_raster,
         p_shape_poly,
         p_shape_raster,
         p_shape)
      gc()

      polygon_forest %>%
        slice(1000) %>%
        sf::st_centroid() %>%
        sf::st_buffer(1000) -> clip_area

      ggplot() +
        # base_map(clip_area %>% st_transform(4326) %>% st_bbox, basemap = "google-satellite", increase_zoom = 3, nolabels = T) +
        geom_sf(data = polygon_forest %>% st_filter(clip_area) %>% rmapshaper::ms_clip(clip_area) %>% sf::st_transform(4326), color = "#D81B60", fill = NA, lwd = 1) +
        geom_sf(data = segmentation_smoothed %>% st_filter(clip_area) %>% rmapshaper::ms_clip(clip_area) %>% sf::st_transform(4326), color = "#1E88E5", fill = NA, lwd = 1) +
        theme_bw() -> p_map

      ggsave(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/maps/map", output_name, ".jpg"),
             plot = p_map,
             width = 10,
             height = 10,
             dpi = 300)

      graphics.off()
      rm(p_map,
         clip_area)
      gc()

      tibble(forest = param$forest,
             model_name = param$model_name,
             thresh = param$thresh,
             spec = param$spec,
             spat = param$spat,

             overlap_area = overlapping::overlap(list(segmentation$area, polygon_forest$area))$OV,
             overlap_area_raster = overlap_area_raster_OV,

             overlap_shape = overlapping::overlap(list(segmentation$shape_index, polygon_forest$shape_index))$OV,
             overlap_shape_raster = overlap_shape_raster_OV,

             overlap_shape_smoothed = overlapping::overlap(list(segmentation_smoothed$shape_index, polygon_forest$shape_index))$OV,
             overlap_shape_smoothed_raster = overlap_shape_smoothed_raster_OV,

             mean_area_forest = mean(polygon_forest$area),
             mean_area_segmentation = mean(segmentation$area),
             mean_area_sspe = sspe(mean_area_forest, mean_area_segmentation),

             sd_area_forest = sd(polygon_forest$area),
             sd_area_segmentation = sd(segmentation$area),
             sd_area_sspe = sspe(sd_area_forest, sd_area_segmentation),

             n_forest = nrow(polygon_forest),
             n_segmentation = nrow(segmentation),
             n_sspe = sspe(n_forest, n_segmentation),

             dens_forest = n_forest/(as.numeric(sum(st_area(polygon_forest)))/1000000),
             dens_segmentation = n_segmentation/(as.numeric(sum(st_area(segmentation)))/1000000),
             dens_sspe = sspe(dens_forest, dens_segmentation)) %>%

        bind_cols(quantile(polygon_forest$area,
                           probs = seq(0, 1, 0.1)) %>%
                    as.list() %>%
                    tibble::as_tibble_row(.name_repair = "minimal") %>%
                    rename_with(~ paste0("quantile_forest_", seq(0, 100, 10))),
                  quantile(segmentation$area,
                           probs = seq(0, 1, 0.1)) %>%
                    as.list() %>%
                    tibble::as_tibble_row(.name_repair = "minimal") %>%
                    rename_with(~ paste0("quantile_segmentation_", seq(0, 100, 10)))) -> results

      results %>%
        saveRDS(paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/", param$forest, "/", param$model_name, "/tables/results", output_name, ".rds"))

      terra::tmpFiles(remove = TRUE)
      rm(metrics,
         polygon_forest,
         segmentation,
         chunk_size,
         segmentation_chunks,
         segmentation_smoothed,
         overlap_area_raster_OV,
         overlap_shape_raster_OV,
         overlap_shape_smoothed_raster_OV)
      gc()

      return(NULL)

    }, error = function(e) {

      cat(c("Error at index ", p, " : ", conditionMessage(e), "----------"),
          file = paste0("D:/00_Ontario_eFRI/segmentation/grm_parametrization/log_error/log_error_", param$forest, "_", param$model_name, "_thresh", param$thresh, "_spec", param$spec*100, "_spat", param$spat*100, ".txt"),
          append = TRUE,
          sep = "\n")

      return(NULL)
    })

  })

}

list.files("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results", full.names = TRUE) %>%
  map(list.files, full.names = TRUE) %>%
  unlist() %>%
  map(list.files, full.names = TRUE, pattern = "tables") %>%
  unlist() %>%
  list.files(full.names = TRUE, pattern = "\\.rds$") %T>%
  {files_names <<- .} %>%
  map_dfr(readRDS) -> results

# Forgot to save the model name in th first round, but its now corrected at line 721, so I need to add it manually for the first results
files_names %>%
  str_extract("model_[A-Z]") -> model_names

results %>%
  mutate(model_name = model_names) %>%
  select(forest, model_name, everything()) -> results

results %>%
  saveRDS("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/results.rds")





# 🟡 Best parameters and figures ----
# Parameters
dens_threshold <- 0.1
weight_area <- 1/2
weight_shape <- 1/2

# Read data
readRDS("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/results.rds") %>%
  mutate(area_score = overlap_area_raster,
         shape_score = overlap_shape_smoothed_raster) %>%
  mutate(final_score = (weight_area*(1-area_score)) + (weight_shape*(1-shape_score))) %>%
  mutate_at(c("thresh", "spec", "spat"), as.numeric) -> results

# Compute for a specific forest and model and produce figures
selected_forest <- "RMF"
selected_forest <- "OVF"
selected_model <- "model_B"

results %>%
  filter(forest == selected_forest) %>%
  filter(model_name == selected_model) %>%
  filter(dens_sspe >= -dens_threshold & dens_sspe <= dens_threshold) %>%
  mutate(parameters = paste0("thresh", thresh, "_spec", spec*100, "_spat", spat*100)) -> results_temp

#plot_ly(data = results_temp %>% filter(spec >= 50),
plot_ly(data = results_temp,
        x = ~area_score, y = ~shape_score, z = ~dens_sspe,
        color = ~final_score,
        colors = c("#117733", "#882235"),
        #size = ~thresh,
        #size = ~spat,
        size = ~spec,
        sizes = c(1, 1000),
        type = "scatter3d",
        mode = "markers",
        text = ~parameters,
        hovertemplate = paste("Overlap Area : %{x}<br>",
                              "Overlap shape : %{y}<br>",
                              "Symmetric Signed Percentage Error (SSPE) : %{z}<br>",
                              "Parameters : %{text}<extra></extra>"))

plot_ly(data = results_temp,
        x = ~thresh, y = ~spec, z = ~spat,
        color = ~ dens_sspe,
        colors = c("#332288", "#117733", "#882235"),
        size = ~ dens_sspe,
        sizes = c(1, 1000),
        type = "scatter3d",
        mode = "markers",
        hovertemplate = paste("Threshold : %{x}<br>",
                              "Spec : %{y}<br>",
                              "Spat : %{z}<br>",
                              "Symmetric Signed Percentage Error (SSPE) : %{marker.color}<extra></extra>"))

plot_ly(data = results_temp,
        x = ~thresh, y = ~spec, z = ~spat,
        color = ~area_score,
        colors = c("#882235", "#117733"),
        size = ~area_score,
        sizes = c(1, 1000),
        type = "scatter3d",
        mode = "markers",
        hovertemplate = paste("Threshold : %{x}<br>",
                              "Spec : %{y}<br>",
                              "Spat : %{z}<br>",
                              "Overlap Area : %{marker.color}<extra></extra>"))

plot_ly(data = results_temp,
        x = ~thresh, y = ~spec, z = ~spat,
        color = ~shape_score,
        colors = c("#882235", "#117733"),
        size = ~shape_score,
        sizes = c(1, 1000),
        type = "scatter3d",
        mode = "markers",
        hovertemplate = paste("Threshold : %{x}<br>",
                              "Spec : %{y}<br>",
                              "Spat : %{z}<br>",
                              "Overlap Shape : %{marker.color}<extra></extra>"))

plot_ly(data = results_temp,
        x = ~thresh, y = ~spec, z = ~spat,
        color = ~final_score,
        colors = c("#117733", "#882235"),
        size = ~final_score,
        sizes = c(1000, 1),
        type = "scatter3d",
        mode = "markers",
        hovertemplate = paste("Threshold : %{x}<br>",
                              "Spec : %{y}<br>",
                              "Spat : %{z}<br>",
                              "Final score : %{marker.color}<extra></extra>"))

results_temp %>%
  arrange(final_score) %>%
  slice(1) %>%
  dplyr::select(forest, model_name, thresh, spec, spat, overlap_area_raster, overlap_shape_smoothed_raster, n_forest, n_segmentation, dens_sspe, final_score) -> best_parameters_temp

# Compute for each forest and model
results %>%
  group_by(forest, model_name) %>%
  reframe() -> forest_model_combinations

map_dfr(seq_len(nrow(forest_model_combinations)),
        function(x){

          results %>%
            filter(forest == forest_model_combinations$forest[x]) %>%
            filter(model_name == forest_model_combinations$model_name[x]) %>%
            filter(dens_sspe >= -dens_threshold & dens_sspe <= dens_threshold) %>%
            # Median values
            # slice_min(abs(spec - median(spec))) %>%
            # slice_min(abs(spat - median(spat))) %>%
            # slice_min(abs(thresh - median(thresh))) %>%
            # Min final score
            arrange(final_score) %>%
            slice(1) %>%
            # Final selection
            dplyr::select(forest, model_name, thresh, spec, spat, overlap_area_raster, overlap_shape_smoothed_raster, n_forest, n_segmentation, dens_sspe, final_score)

        }) %>%
  mutate(type = "by_forest_model") -> best_parameters_by_forest_model

# Compute for each forest
best_parameters_by_forest <- results %>%
  group_by(forest, thresh, spec, spat) %>%
  reframe(area_score = mean(area_score),
          shape_score = mean(shape_score),
          dens_sspe = mean(dens_sspe)) %>%
  mutate(final_score = (weight_area*(1-area_score)) + (weight_shape*(1-shape_score))) %>%
  filter(dens_sspe >= -dens_threshold & dens_sspe <= dens_threshold) %>%
  group_by(forest) %>%
  # # Median values
  # slice_min(abs(spec - median(spec))) %>%
  # slice_min(abs(spat - median(spat))) %>%
  # slice_min(abs(thresh - median(thresh)))
  # Min final score
  arrange(final_score) %>%
  slice(1)

map_dfr(unique(results$forest),
        function(f){

          best_parameters_by_forest %>%
            filter(forest == f) -> best_parameters_by_forest_temp

          results %>%
            filter(forest == f) %>%
            filter(thresh == best_parameters_by_forest_temp$thresh &
                     spec == best_parameters_by_forest_temp$spec &
                     spat == best_parameters_by_forest_temp$spat) %>%
            dplyr::select(forest, model_name, thresh, spec, spat, overlap_area_raster, overlap_shape_smoothed_raster, n_forest, n_segmentation, dens_sspe, final_score)

        }) %>%
  mutate(type = "by_forest") -> best_parameters_by_forest

# Compute general results
best_parameters <- results %>%
  group_by(thresh, spec, spat) %>%
  reframe(area_score = mean(area_score),
          shape_score = mean(shape_score),
          dens_sspe = mean(dens_sspe)) %>%
  mutate(final_score = (weight_area*(1-area_score)) + (weight_shape*(1-shape_score))) %>%
  filter(dens_sspe >= -dens_threshold & dens_sspe <= dens_threshold) %>%
  # # Median values
  # slice_min(abs(spec - median(spec))) %>%
  # slice_min(abs(spat - median(spat))) %>%
  # slice_min(abs(thresh - median(thresh)))
  # # Min final score
  arrange(final_score) %>%
  slice(1)

results %>%
  filter(thresh == best_parameters$thresh &
           spec == best_parameters$spec &
           spat == best_parameters$spat) %>%
  dplyr::select(forest, model_name, thresh, spec, spat, overlap_area_raster, overlap_shape_smoothed_raster, n_forest, n_segmentation, dens_sspe, final_score) %>%
  mutate(type = "all") -> best_parameters

# Merge all results
bind_rows(best_parameters_by_forest_model,
          best_parameters_by_forest,
          best_parameters) -> best_parameters_all

best_parameters_all %>%
  saveRDS("D:/00_Ontario_eFRI/segmentation/grm_parametrization/results/best_parameters_all.rds")
