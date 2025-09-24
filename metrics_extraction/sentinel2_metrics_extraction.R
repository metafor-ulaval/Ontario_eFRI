# 🟡 Library ----
library(tidyverse)
library(terra)





# 🟡 Setwd and parameters ----
sector <- "OVF"
setwd(paste0("D:/00_Ontario_eFRI/", sector))





# 🟡 Read data ----
# rast("OV_20m.tif") -> sentinel2_2
# 
# rast("OV_10m.tif") %>% 
#   terra::resample(sentinel2_2, method = "average") -> sentinel2_1
# 
# c(sentinel2_1,
#   sentinel2_2) -> sentinel2
# 
# sentinel2 %>% 
#   names %>% 
#   gsub("_median", ., replacement = "") -> names(sentinel2)
# 
# names(sentinel2)
# 
# sentinel2 %>% 
#   writeRaster("sentinel_all_bands_20m.tif")
              
rast("./metrics/sentinel2/RAW/sentinel_all_bands_20m.tif") -> sentinel_2_20m





# 🟡 Export each band ----
map(names(sentinel_2_20m), ~{
  sentinel_2_20m[[.x]] %>% 
    writeRaster(paste0("./metrics/sentinel2/", .x, ".tif"))
})





# 🟡 Normalized differences veg. i. (NDVI) ----
rast("./metrics/sentinel2/B4.tif") -> B4
rast("./metrics/sentinel2/B8.tif") -> B8
(B8 - B4)/(B8 + B4) -> NDVI
names(NDVI) <- "NDVI"
NDVI %>% 
  writeRaster("./metrics/sentinel2/NDVI.tif")





# 🟡 Green normalized differences veg. i. (GNDVI) ----
rast("./metrics/sentinel2/B3.tif") -> B3
rast("./metrics/sentinel2/B8.tif") -> B8
(B8 - B3)/(B8 + B3) -> GNDVI
names(GNDVI) <- "GNDVI"
GNDVI %>% 
  writeRaster("./metrics/sentinel2/GNDVI.tif")





# 🟡 Normalized differences veg. i. red-egde (NDVIre) ----
rast("./metrics/sentinel2/B5.tif") -> B5
rast("./metrics/sentinel2/B8.tif") -> B8
(B8 - B5)/(B8 + B5) -> NDVIre
names(NDVIre) <- "NDVIre"
NDVIre %>% 
  writeRaster("./metrics/sentinel2/NDVIre.tif")





# 🟡 Ratio veg. i. (RVI) ----
rast("./metrics/sentinel2/B4.tif") -> B4
rast("./metrics/sentinel2/B8.tif") -> B8
B8/B4 -> RVI
names(RVI) <- "RVI"
RVI %>% 
  writeRaster("./metrics/sentinel2/RVI.tif")





# 🟡 Difference veg. i. (DVI) ----
rast("./metrics/sentinel2/B4.tif") -> B4
rast("./metrics/sentinel2/B8.tif") -> B8
B8-B4 -> DVI
names(DVI) <- "DVI"
DVI %>% 
  writeRaster("./metrics/sentinel2/DVI.tif")





# 🟡 Soil-adjusted veg. i. (SAVI) ----
rast("./metrics/sentinel2/B4.tif") -> B4
rast("./metrics/sentinel2/B8.tif") -> B8
L <- 0.5
((B8-B4)/(B8+B4+L))*(1+L) -> SAVI
names(SAVI) <- "SAVI"
SAVI %>% 
  writeRaster("./metrics/sentinel2/SAVI.tif")





# 🟡 Enhanced veg. i. (EVI)  ----
rast("./metrics/sentinel2/B2.tif") -> B2
rast("./metrics/sentinel2/B4.tif") -> B4
rast("./metrics/sentinel2/B8.tif") -> B8
2.5*(B8-B4)/(B8+6*B4-7.5*B2+1) -> EVI
names(EVI) <- "EVI"
EVI %>% 
  writeRaster("./metrics/sentinel2/EVI.tif")





# 🟡 Green chlorophyll i. (GCI) ----
rast("./metrics/sentinel2/B3.tif") -> B3
rast("./metrics/sentinel2/B8.tif") -> B8
(B8/B3)-1 -> GCI
names(GCI) <- "GCI"
GCI %>% 
  writeRaster("./metrics/sentinel2/GCI.tif")





# 🟡 Normalized difference moisture i. (NDMI)  ----
rast("./metrics/sentinel2/B8.tif") -> B8
rast("./metrics/sentinel2/B11.tif") -> B11
(B8-B11)/(B8+B11) -> NDMI
names(NDMI) <- "NDMI"
NDMI %>% 
  writeRaster("./metrics/sentinel2/NDMI.tif")





# 🟡 Normalized difference water i. (NDWI)  ----
rast("./metrics/sentinel2/B8.tif") -> B8
rast("./metrics/sentinel2/B12.tif") -> B12
(B8-B12)/(B8+B12) -> NDWI
names(NDWI) <- "NDWI"
NDWI %>% 
  writeRaster("./metrics/sentinel2/NDWI.tif")