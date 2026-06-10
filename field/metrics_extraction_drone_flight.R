# 🟡 Library ----
library(tidyverse)
library(terra)
library(sf)
library(lasR)
library(lidR)
library(lidRmetrics)
library(future)
library(RSAGA)
library(whitebox)
RSAGA::rsaga.env(path = "C:/Logiciels/saga-8.3.0_x64/saga-8.3.0_x64") -> env_saga





# 🟡 Parameters and functions ----
source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/functions_archives/other/split_metrics.R")
plan(multisession, workers = 5L)
set_lidr_threads(0)

# sub_sector <- "S1"





# 🟡 Buffered Area Based Approach ----
for(sub_sector in paste0("S", 1:10)){

  # Message
  cat(paste0("Computing LiDAR metrics for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/metrics"))

  # Read catalog
  readLAScatalog(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/nuages/02_norm_decimated")) -> ctg_normalized

  # Read steps
  source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/Ontario_eFRI/metrics_extraction/lasR_metrics_list.R")
  source("C:/Users/FRLES121/OneDrive - Université Laval/Documents/Ontario_eFRI/Scripts/Interne/Ontario_eFRI/metrics_extraction/lidRmetrics_metrics_list.R")

  read <- reader(filter = "-keep_random_fraction 0.01")

  start <- Sys.time()
  # Build lasR pipe
  ls(pattern = "lasR") %>%
    c("read", .) %>%
    paste0(collapse = " + ") %>%
    parse(text = .) %>%
    eval() -> pipeline_1

  # Apply lasR pipe
  lasR::exec(pipeline_1,
             on = ctg_normalized,
             progress = TRUE,
             ncores = nested(ncores = ceiling((ncores()-4L)/8L), ncores2 = 8L))

  # Build lidRmetrics pipe
  ls(pattern = "lidRmetrics") %>%
    c("read", .) %>%
    paste0(collapse = " + ") %>%
    parse(text = .) %>%
    eval() -> pipeline_2

  # Apply lidRmetrics pipe
  lasR::exec(pipeline_2,
             on = ctg_normalized,
             progress = TRUE,
             ncores = nested(ncores = ceiling((ncores()-4L)/8L), ncores2 = 8L))
  end <- Sys.time()
  duration <- end - start
  duration
}





# 🟡 Split metrics ----
for(sub_sector in paste0("S", 1:10)){
  # Message
  cat(paste0("Split metrics for sector : ", sub_sector, "\n"))

  # Computing
  split_metrics(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/metrics"))
}





# 🟡 DEM ----
for(sub_sector in paste0("S", 1:10)){

  # Message
  cat(paste0("Computing DEM for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/metrics"))

  # Read catalog
  readLAScatalog(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/nuages/03_decimated")) -> ctg

  # DEM 1 m
  read <- reader()

  dem_1m <- lasR::dtm(1,
                      ofile = "./dem_temp_1m.tif")

  # Apply pipe
  lasR::exec(read + dem_1m,
             on = ctg,
             progress = TRUE,
             ncores = nested(ncores = ceiling((ncores()-4)/4L), ncores2 = 4L))

  # Rewrite dem
  rast("./dem_temp_1m.tif") %>%
    writeRaster("./dem_1m.tif")

  file.remove("./dem_temp_1m.tif")
}





# 🟡 DEM buffered ----
for(sub_sector in paste0("S", 1:10)){

  # Message
  cat(paste0("Computing dem buffered for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("D:/dossier_remise_fin_contrat_PY/ontario/metriques/", sub_sector, "/metrics"))

  rast("./dem_1m.tif") -> dem
  dem %>%
    terra::focal(w = 21,
                 fun = "mean",
                 na.policy = "omit",
                 na.rm = TRUE,
                 expand = TRUE) -> dem_buffered

  names(dem_buffered) <- "dem_buffered"

  dem_buffered %>%
    writeRaster("./dem_buffered.tif")

}







# 🟡 Slope ----
for(sub_sector in paste0("S", 1:10)){

  # Message
  cat(paste0("Computing slope for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("D:/dossier_remise_fin_contrat_PY/ontario/metriques/", sub_sector, "/metrics"))

  rast("./dem_1m.tif") -> dem
  dem %>%
    terrain(v = "slope",
            unit = "radians") %>%
    tan() %>%
    {.*100} -> slope

  slope %>%
    writeRaster("./slope.tif")

}





# 🟡 Slope buffered ----
for(sub_sector in paste0("S", 1:10)){

  # Message
  cat(paste0("Computing slope buffered for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("D:/dossier_remise_fin_contrat_PY/ontario/metriques/", sub_sector, "/metrics"))

  rast("./slope.tif") -> slope
  slope %>%
    terra::focal(w = 21,
                 fun = "mean",
                 na.policy = "omit",
                 na.rm = TRUE,
                 expand = TRUE) -> slope_buffered

  names(slope_buffered) <- "slope_buffered"

  slope_buffered %>%
    writeRaster("./slope_buffered.tif")

}







# 🟡 Merge DEM ----
for(sub_sector in paste0("S", 1:10)){

  # Message
  cat(paste0("Merging DEM for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/metrics"))

  rast("./dem_1m.tif") -> dem_1m

  rast("./dem_large_1m.tif") -> dem_large_1m

  dem_large_1m %>%
    crop(dem_1m) %>%
    mask(dem_1m) -> dem_large_1m_crop

  global((dem_large_1m_crop - dem_1m),
         "mean",
         na.rm = TRUE) %>%
    pull(mean) -> mean

  dem_1m + mean -> dem_1m_corrected

  project(dem_1m_corrected, dem_large_1m) -> dem_1m_corrected

  cover(dem_1m_corrected, dem_large_1m) %>%
    writeRaster("./dem_large_merged_1m.tif")

}






# 🟡 SAGAWI ----
# sub_sector <- "S1"
# dem <- "dem_1m.tif"

for(sub_sector in paste0("S", 6)){

  # Message
  cat(paste0("Computing SAGAWI for sector : ", sub_sector, "\n"))

  for(dem in c("dem_1m.tif",
               "dem_large_1m.tif",
               "dem_large_merged_1m.tif")){

    cat(paste0("Dem used : ", dem, "\n"))

    # Setwd
    setwd(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/metrics"))

    # SAGAWI
    wbt_feature_preserving_smoothing(dem = dem,
                                     output = "./dem_smoothed_1m.tif")

    wbt_breach_depressions_least_cost(dem = "./dem_smoothed_1m.tif",
                                      output = "./dem_breach_lc_1m.tif",
                                      dist = 50,
                                      flat_increment = 0.0001,
                                      fill = FALSE)

    wbt_breach_depressions(dem = "./dem_breach_lc_1m.tif",
                           output = "./dem_breach_f_1m.tif",
                           flat_increment = 0.0001,
                           fill_pits = TRUE)

    sagawi_temp_output <- gsub("dem", dem, replacement = "sagawi_temp")

    RSAGA::rsaga.wetness.index(in.dem = "./dem_breach_f_1m.tif",
                               out.wetness.index = sagawi_temp_output,
                               suction = 10,
                               area.type = "absolute",
                               slope.type = "local",
                               env = env_saga)

    rast(sagawi_temp_output) -> sagawi_temp

    crs(sagawi_temp) <- paste0("epsg:", crs(rast(dem), describe = TRUE)$code)

    sagawi_output <- gsub("dem", dem, replacement = "sagawi")

    sagawi_temp %>%
      writeRaster(sagawi_output)

    file.remove("./dem_smoothed_1m.tif",
                "./dem_breach_lc_1m.tif",
                "./dem_breach_f_1m.tif",
                sagawi_temp_output)
  }
}
