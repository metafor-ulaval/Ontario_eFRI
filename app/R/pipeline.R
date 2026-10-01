# 🟡🟡 eFRI pipeline 🟡🟡 ----
# This file is sourced in a background R process launched by app.R (callr::r_bg).
# Everything printed with cat() goes to the log file of the segmentation.
# The current step is written in a progress file read by the app every second.





# 🟡🟡 Functions 🟡🟡 ----
# 🌐 Write current step in progress file 🌐
write_progress <- function(progress_file, step, n_steps, message) {
  writeLines(c(step, n_steps, message), progress_file)
  cat(paste0("\n==== Step ", step, "/", n_steps, " : ", message, " ====\n"))
}

# 🌐 Read a raster on the template grid, cropping before projecting 🌐
# When an area is given, the source raster is cropped to the area (with a buffer to keep
# the neighbours used by bilinear interpolation) before being projected on the template.
read_raster_on_template <- function(path,
                                    template,
                                    area = NULL,
                                    method = "bilinear") {
  r <- terra::rast(path)

  if (!is.null(area)) {
    area %>%
      terra::buffer(4 * terra::res(template)[1]) %>%
      terra::project(terra::crs(r)) -> area_source_crs

    r <- terra::crop(r, area_source_crs, snap = "out")
  }

  r <- terra::project(r, template, method = method)

  if (!is.null(area)) {
    r <- terra::mask(r, area)
  }

  return(r)
}

# 🌐 Read a vector layer, only reading features of the area when one is given 🌐
read_vector_in_area <- function(path,
                                crs_target,
                                area = NULL) {
  if (is.null(area)) {
    return(sf::st_transform(sf::st_read(path, quiet = TRUE), crs_target))
  }

  layer_crs <- sf::st_crs(sf::st_layers(path)$crs[[1]])

  area %>%
    sf::st_union() %>%
    sf::st_transform(layer_crs) %>%
    sf::st_as_text() -> area_wkt

  sf::st_read(path, quiet = TRUE, wkt_filter = area_wkt) %>%
    sf::st_transform(crs_target) %>%
    sf::st_filter(area)
}





# 🟡🟡 Pipeline 🟡🟡 ----
run_efri_pipeline <- function(p) {

  suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(sf)
    library(terra)
    library(eFRItools)
  })

  n_steps <- 7

  # 🟢 Memory 🟢
  Sys.setenv(OTB_MAX_RAM_HINT = p$ram_mb)
  Sys.setenv(OTB_MEMORY_AVAILABLE = p$ram_mb)
  Sys.setenv(GDAL_CACHEMAX = p$gdal_cache_mb)

  segmentation_wd <- p$segmentation_wd
  metadata_file <- paste0(segmentation_wd, "/metadata.txt")

  add_metadata <- function(title, value) {
    cat(c(title, value, "----------"),
        file = metadata_file,
        append = TRUE,
        sep = "\n")
  }

  # 🟢 Template grid and area 🟢
  read_metrics(list.files(paste0(p$wd, "/metrics/", p$forest), full.names = TRUE)) -> metrics_infos

  metrics_infos %>%
    filter(name == "dem") %>%
    pull(path) %>%
    rast() -> epsg_rast

  epsg_rast %>%
    st_crs() -> epsg

  if (!is.null(p$extraction_area)) {
    p$extraction_area %>%
      st_transform(epsg) -> extraction_area

    extraction_area_vect <- vect(extraction_area)
    template <- crop(epsg_rast, extraction_area_vect)
  } else {
    extraction_area <- NULL
    extraction_area_vect <- NULL
    template <- epsg_rast
  }

  # 🟢 Best models for imputation 🟢
  list.files(paste0(p$wd, "/analysis/imputation"), pattern = paste0("results_", p$forest), full.names = TRUE) %>%
    map_dfr(function(x){
      read.csv(x) %>%
        arrange(desc(accuracy)) %>%
        slice_max(accuracy, n = 1)
    }) -> imputation_results

  imputation_results %>%
    pull(knn_vars) %>%
    str_remove(",X,Y") %>%
    strsplit(",") %>%
    unlist() %>%
    unique() -> imputation_metrics

  # 🟢 Save metadata 🟢
  start <- Sys.time()
  add_metadata("Start time : ", as.character(start))
  add_metadata("Segmentation metrics : ", p$segmentation_metrics)
  add_metadata("Summary metrics : ", p$summary_metrics)
  add_metadata("EPSG : ", epsg$input)
  add_metadata("Segmentation parameters : ", paste0("thresh = ", p$thresh, " / spec = ", p$spec, " / spat = ", p$spat))
  add_metadata("Segmentation masks : ", paste0(p$masks, collapse = ", "))
  add_metadata("Memory (MB) : ", paste0("OTB = ", p$ram_mb, " / GDAL cache = ", p$gdal_cache_mb))

  if (!is.null(extraction_area)) {
    extraction_area %>%
      st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
               layer = "extraction_area",
               quiet = TRUE)
  }

  # 🟢 Read metrics 🟢
  write_progress(p$progress_file, 1, n_steps, "Read metrics")

  metrics_infos %>%
    filter(name %in% c(p$segmentation_metrics,
                       imputation_metrics,
                       p$summary_metrics,
                       "z_p95", "z_above2", "slope", "sagawi")) -> metrics_infos_selected

  metrics_infos_selected %>%
    pull(path) %>%
    map(read_raster_on_template, template = template, area = extraction_area_vect, method = "bilinear") %>%
    rast() -> metrics

  names(metrics) <- pull(metrics_infos_selected, name)

  # 🟢 Read vector data 🟢
  write_progress(p$progress_file, 2, n_steps, "Read catalog, masks and forest inventory polygons")

  # Catalog
  read_vector_in_area(paste0(p$wd, "/shapefiles/", p$forest, "/ctg.shp"), epsg, extraction_area) -> ctg

  # Masks
  if (!is.null(p$masks)) {
    paste0(p$wd, "/shapefiles/", p$forest, "/", p$masks, ".shp") %>%
      map(function(path){
        if (!is.null(extraction_area_vect)) {
          layer_crs <- crs(vect(path, proxy = TRUE))
          vect(path, filter = project(extraction_area_vect, layer_crs)) %>%
            project(epsg_rast) %>%
            crop(extraction_area_vect) %>%
            mask(extraction_area_vect)
        } else {
          project(vect(path), epsg_rast)
        }
      }) -> masks
  } else {
    masks <- NULL
  }

  # Forest inventory polygons (fri)
  read_vector_in_area(paste0(p$wd, "/shapefiles/", p$forest, "/PolygonForest.shp"), epsg, extraction_area) %>%
    rowid_to_column("id") -> fri_polygons

  # 🟢 Read landcover, forest age and disturbances 🟢
  write_progress(p$progress_file, 3, n_steps, "Read landcover, forest age, forest fire and forest harvest")

  read_raster_on_template(paste0(p$wd, "/metrics/", p$forest, "/other/landcover.tif"),
                          template, extraction_area_vect, method = "near") -> landcover

  landcover_codes <- c(`1` = "Clear_Open_Water",
                       `2` = "Turbid_Water",
                       `3` = "Shoreline",
                       `4` = "Mudflats",
                       `5` = "Marsh",
                       `6` = "Swamp",
                       `7` = "Fen",
                       `8` = "Bog",
                       `10` = "Heath",
                       `11` = "Sparse_Treed",
                       `12` = "Treed_Upland",
                       `13` = "Deciduous_Treed",
                       `14` = "Mixed_Treed",
                       `15` = "Coniferous_Treed",
                       `16` = "Plantations_Treed_Cultivated",
                       `17` = "Hedge_Rows",
                       `18` = "Disturbance",
                       `19` = "Open_Cliff_Talus",
                       `20` = "Alvar",
                       `21` = "Sand_Barren_Dune",
                       `22` = "Open_Tallgrass_Prairie",
                       `23` = "Tallgrass_Savannah",
                       `24` = "Tallgrass_Woodland",
                       `25` = "Sand_Gravel_Mine_Tailings_Extraction",
                       `26` = "Bedrock",
                       `27` = "Communit_Infrastructure",
                       `28` = "Agriculture_Undifferentiated_Rural_Land_Use",
                       `157` = "Other",
                       `247` = "Cloud_Shadow")

  # Attach levels (lookup table)
  levels(landcover) <- data.frame(value = as.integer(names(landcover_codes)),
                                  class = unname(landcover_codes))

  read_raster_on_template(paste0(p$wd, "/metrics/", p$forest, "/other/forest_age_2019.tif"),
                          template, extraction_area_vect, method = "near") -> forest_age_2019

  read_raster_on_template(paste0(p$wd, "/metrics/", p$forest, "/other/forest_fire_1985_2020.tif"),
                          template, extraction_area_vect, method = "near") -> forest_fire_1985_2020

  read_raster_on_template(paste0(p$wd, "/metrics/", p$forest, "/other/forest_harvest_1985_2020.tif"),
                          template, extraction_area_vect, method = "near") -> forest_harvest_1985_2020

  # 🟢 Segmentation 🟢
  write_progress(p$progress_file, 4, n_steps, "Perform segmentation")

  eFRI_segmentation(metrics = metrics[[p$segmentation_metrics]],
                    masks = masks,
                    thresh = p$thresh,
                    spec = p$spec,
                    spat = p$spat,
                    method = "bs",
                    clean_nodata = TRUE,
                    output_path = segmentation_wd,
                    output_name = "segmentation",
                    otb_dir = p$otb_dir) -> segmentation

  segmentation %>%
    st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
             layer = "segmentation",
             quiet = TRUE)

  # 🟢 Build attribute table 🟢
  write_progress(p$progress_file, 5, n_steps, "Build attribute table of segmented polygons")

  eFRI_attribute_table(segmentation = segmentation,
                       metrics = metrics,
                       summary_metrics = unique(c(p$summary_metrics, "z_p95", "z_above2", "slope", "sagawi")),
                       landcover = landcover,
                       forest_fire = forest_fire_1985_2020,
                       forest_harvest = forest_harvest_1985_2020,
                       forest_age = forest_age_2019) -> segmentation_data

  segmentation_data %>%
    dplyr::select(-id_seg, -id) %>%
    rename(HEIGHT = Z_P95,
           CANOPY_COVER = Z_ABOVE2,
           MOISTURE = SAGAWI) -> segmentation_data

  segmentation_data %>%
    st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
             layer = "segmentation_data",
             quiet = TRUE)

  # 🟢 Imputation 🟢
  write_progress(p$progress_file, 6, n_steps, "Perform imputation")

  eFRI_imputation(segmentation = segmentation_data,
                  forest_polygon = fri_polygons,
                  metrics = metrics,
                  landcover = landcover,
                  forest_fire = forest_fire_1985_2020,
                  forest_harvest = forest_harvest_1985_2020,
                  ctg = ctg,
                  lidar_year_field = "Fl_Cr_Y",
                  forest_year_field = "YRUPD",
                  forest_composition_field = "SPCOMP",
                  forest_type_field = "POLYTYPE",
                  target_var = imputation_results$target_var,
                  knn_var = imputation_results$knn_vars) -> segmentation_data_imputed

  # 🟢 Write outputs and metadata 🟢
  write_progress(p$progress_file, 7, n_steps, "Write outputs")

  segmentation_data_imputed %>%
    st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
             layer = "segmentation_data_imputed",
             quiet = TRUE)

  add_metadata("Number of forest polygon created : ", nrow(segmentation_data))

  end <- Sys.time()
  add_metadata("End time : ", as.character(end))

  elapsed_time <- sub("Time difference of ", "", capture.output(difftime(end, start)))
  cat(c("Elapsed time : ", elapsed_time),
      file = metadata_file,
      append = TRUE,
      sep = "\n")

  return(list(n_polygons = nrow(segmentation_data),
              elapsed_time = elapsed_time))
}
