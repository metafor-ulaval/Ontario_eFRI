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





# 🟡 Dem and SAGAWI ----
for(sub_sector in paste0("S", 2:10)){

  # Message
  cat(paste0("Computing DEM and SAGAWI for sector : ", sub_sector, "\n"))

  # Setwd
  setwd(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/metrics"))

  # Read catalog
  readLAScatalog(paste0("F:/Ontario_eFRI/02_donnees_traitees/", sub_sector, "/nuages/03_decimated")) -> ctg

  # DEM 1 m
  read <- reader()

  dem_1m <- lasR::dtm(1,
                      ofile = "./dem_1m.tif")

  # Apply pipe
  lasR::exec(read + dem_1m,
             on = ctg,
             progress = TRUE,
             ncores = nested(ncores = ceiling((ncores()-4)/4L), ncores2 = 4L))

  # Rewrite dem
  rast("./dem_1m.tif") %>%
    writeRaster("./dem_1m_rewrite.tif")

  # SAGAWI 1 m
  wbt_feature_preserving_smoothing(dem = "./dem_1m_rewrite.tif",
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

  RSAGA::rsaga.wetness.index(in.dem = "./dem_breach_f_1m.tif",
                             out.wetness.index = "./SAGAWI_1m.tif",
                             suction = 10,
                             area.type = "absolute",
                             slope.type = "local",
                             env = env_saga)
}
