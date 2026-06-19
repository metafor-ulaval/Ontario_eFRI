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
RSAGA::rsaga.env(path = "C:/Logiciels/saga-9.12.5_msw") -> env_saga





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
                      ofile = "./dem.tif")

  # Apply pipe
  lasR::exec(read + dem_1m,
             on = ctg,
             progress = TRUE,
             ncores = nested(ncores = ceiling((ncores()-4)/4L), ncores2 = 4L))

}





# 🟡 Other topo ----
x <- "D:/00_Ontario_eFRI/random_forest/metriques/S6"

list.files("D:/00_Ontario_eFRI/random_forest/metriques", full.names = T) %>%
  map(function(x){

    # Dem
    rast(paste0(x, "/metrics/temp/dem.tif")) %>% terra::focal(w = 21, fun = "mean") -> dem
    names(dem) <- "dem"
    dem %>% writeRaster(paste0(x, "/metrics/dem.tif"))

    # Slope
    rast(paste0(x, "/metrics/temp/dem.tif")) %>% terrain(v = "slope", unit = "degrees") -> slope

    names(slope) <- "slope"
    slope %>% writeRaster(paste0(x, "/metrics/temp/slope.tif"))

    slope %>% terra::focal(w = 21, fun = "mean") -> slope
    names(slope) <- "slope"
    slope %>% writeRaster(paste0(x, "/metrics/slope.tif"))

    # Aspect
    rast(paste0(x, "/metrics/temp/dem.tif")) %>% terrain(v = "aspect", unit = "degrees") -> aspect

    names(aspect) <- "aspect"
    aspect %>% writeRaster(paste0(x, "/metrics/temp/aspect.tif"))

    sin(aspect*pi/180) -> eastness
    names(eastness) <- "eastness"
    eastness %>% writeRaster(paste0(x, "/metrics/temp/eastness.tif"))

    eastness %>% terra::focal(w = 21, fun = "mean") -> eastness
    names(eastness) <- "eastness"
    eastness %>% writeRaster(paste0(x, "/metrics/eastness.tif"))

    cos(aspect*pi/180) -> northness
    names(northness) <- "northness"
    northness %>% writeRaster(paste0(x, "/metrics/temp/northness.tif"))

    northness %>% terra::focal(w = 21, fun = "mean") -> northness
    names(northness) <- "northness"
    northness %>% writeRaster(paste0(x, "/metrics/northness.tif"))

    # Breach depression
    whitebox::wbt_breach_depressions_least_cost(dem = paste0(x, "/metrics/temp/dem.tif"),
                                                output = paste0(x, "/metrics/temp/dem_breach_lc.tif"),
                                                dist = 50,
                                                flat_increment = 0.0001,
                                                fill = FALSE)

    whitebox::wbt_breach_depressions(dem = paste0(x, "/metrics/temp/dem_breach_lc.tif"),
                                     output =  paste0(x, "/metrics/temp/dem_breach_f.tif"),
                                     flat_increment = 0.0001,
                                     fill_pits = TRUE)

    # TWI
    whitebox::wbt_fd8_flow_accumulation(dem = paste0(x, "/metrics/temp/dem_breach_f.tif"),
                                        output = paste0(x, "/metrics/temp/facc_d8.tif"))

    whitebox::wbt_slope(dem = paste0(x, "/metrics/temp/dem_breach_f.tif"),
                        output = paste0(x, "/metrics/temp/slope_wbt.tif"))

    whitebox::wbt_wetness_index(sca = paste0(x, "/metrics/temp/facc_d8.tif"),
                                slope = paste0(x, "/metrics/temp/slope_wbt.tif"),
                                output = paste0(x, "/metrics/temp/twi.tif"))

    rast(paste0(x, "/metrics/temp/twi.tif")) %>% terra::focal(w = 21, fun = "mean") -> twi
    names(twi) <- "twi"
    twi %>% writeRaster(paste0(x, "/metrics/twi.tif"))

    # SAGAWI
    RSAGA::rsaga.wetness.index(in.dem = paste0(x, "/metrics/temp/dem_breach_f.tif"),
                               out.wetness.index = paste0(x, "/metrics/temp/sagawi.tif"),
                               suction = 10,
                               area.type = "absolute",
                               slope.type = "local",
                               env = env_saga)

    rast(paste0(x, "/metrics/temp/sagawi.tif")) %>% terra::focal(w = 21, fun = "mean") -> sagawi
    names(sagawi) <- "sagawi"
    sagawi %>% writeRaster(paste0(x, "/metrics/sagawi.tif"))

    # Relative topographic position
    weiss_topographic_position_index(dem = rast(paste0(x, "/metrics/temp/dem_breach_f.tif")),
                                     inner_radius = 25,
                                     outer_radius = 50,
                                     round = FALSE) -> tpi

    names(tpi) <- "tpi"
    tpi %>% writeRaster(paste0(x, "/metrics/temp/tpi.tif"))

    tpi %>% terra::focal(w = 21, fun = "mean") -> tpi
    names(tpi) <- "tpi"
    tpi %>% writeRaster(paste0(x, "/metrics/tpi.tif"))

  })
