# 🟡 Library ----
library(tidyverse)
library(whitebox)
library(terra)
library(sf)
library(lasR)
library(lidR)
library(RSAGA)
library(lidRmetrics)





# 🟡 Setwd, parameters and functions ----
sector <- "OVF"
setwd(paste0("D:/00_Ontario_eFRI/", sector))





# 🟡 Metrics extraction ----
# Landcover
rast("D:/00_Ontario_eFRI/data/landcover/OLCC_V2_TIFF.tif") %>% 
  terra::project(rast("./metrics/lidar/dem.tif"), method = "near") %>% 
  mask(rast("./metrics/lidar/dem.tif")) -> landcover

landcover %>% 
  writeRaster("./metrics/other/landcover.tif")


# Forest_Fire_1985-2020
rast("D:/00_Ontario_eFRI/data/CA_Forest_Fire_1985-2020/CA_Forest_Fire_1985-2020.tif") %>% 
  terra::project(rast("./metrics/lidar/dem.tif"), method = "near") %>% 
  mask(rast("./metrics/lidar/dem.tif")) -> forest_fire_1985_2020

forest_fire_1985_2020 %>% 
  writeRaster("./metrics/other/forest_fire_1985_2020.tif")


# Forest_Harvest_1985-2020
rast("D:/00_Ontario_eFRI/data/CA_Forest_Harvest_1985-2020/CA_Forest_Harvest_1985-2020.tif") %>% 
  terra::project(rast("./metrics/lidar/dem.tif"), method = "near") %>% 
  mask(rast("./metrics/lidar/dem.tif")) -> forest_harvest_1985_2020

forest_harvest_1985_2020 %>%
  writeRaster("./metrics/other/forest_harvest_1985_2020.tif")


# Forest_Age_2019
rast("D:/00_Ontario_eFRI/data/CA_forest_age_2019/CA_forest_age_2019.tif") %>% 
  terra::project(rast("./metrics/lidar/dem.tif"), method = "near") %>% 
  mask(rast("./metrics/lidar/dem.tif")) -> forest_age_2019

forest_age_2019 %>% 
  writeRaster("./metrics/other/forest_age_2019.tif")