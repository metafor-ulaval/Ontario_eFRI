# 🟡 Library ----
library(tidyverse)
library(whitebox)
library(terra)
library(sf)
library(lasR)
library(lidR)
library(RSAGA)
library(lidRmetrics)
library(future)





# 🟡 Setwd, parameters and functions ----
sector <- "RMF"
setwd(paste0("D:/00_Ontario_eFRI/", sector))

RSAGA::rsaga.env(path = "D:/00_Ontario_eFRI/logiciels/saga-9.6.2_x64") -> env_saga

# Cannot use parallel computing
# fractional.cover <- function(z,
#                              lower,
#                              upper){
#   return(sum(z > lower & z < upper)/sum(z < upper))
# }

# dem




# # 🟡 Create calatog ----
# Create index (.lax)
# list.files(path = "./laz",
#            full.names = T) %>%
#   map(rlas::writelax)
#
# # Create catalog
# readLAScatalog("./laz") %>%
#   saveRDS("./ctg/ctg.rds")
#
# # Create shapefile
# readRDS("./ctg/ctg.rds") %>%
#   st_as_sf() %>%
#   st_write("./ctg/ctg.shp")
#
#
#
#
#
# 🟡 Read calatog ----
readRDS("./ctg/ctg.rds") -> ctg

las_check(ctg)

ctg@data$CRS %>%
  table

ctg@data$Number.of.variable.length.records %>%
  table

opt_chunk_buffer(ctg) <- 50

plot(ctg, chunk_pattern = T)





# # 🟡 Test petit ctg ----
# st_read("./ctg/ctg.shp") %>%
#   slice(2000) %>%
#   st_buffer(2000) %>%
#   st_geometry() %>%
#   st_filter(st_read("./ctg/ctg.shp") ,.) %>%
#   pull(filenam) %>%
#   readLAScatalog() -> ctg
#
#
#
#
#
#
# # 🟡 Tri calatog OVF ----
# st_read("./ctg/ctg.shp") %>%
#   filter(Nmb____ == 1) %>%
#   pull(filenam) %>%
#   readLAScatalog() %>%
#   saveRDS("./ctg/ctg2.rds")
#
#
#
#
#
# 🟡 Metrics from laz data ----
# Create pipe steps
read <- reader()

normalize <- lasR::normalize()

dem <- lasR::dtm(20,
                 ofile = "./metrics/lidar/dem.tif")

z_p5 <- lasR::rasterize(20,
                        operators = "z_p5",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p5.tif")

z_p10 <- lasR::rasterize(20,
                        operators = "z_p10",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p10.tif")

z_p20 <- lasR::rasterize(20,
                        operators = "z_p20",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p20.tif")

z_p30 <- lasR::rasterize(20,
                        operators = "z_p30",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p30.tif")

z_p40 <- lasR::rasterize(20,
                        operators = "z_p40",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p40.tif")

z_p50 <- lasR::rasterize(20,
                        operators = "z_p50",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p50.tif")

z_p60 <- lasR::rasterize(20,
                        operators = "z_p60",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p60.tif")

z_p70 <- lasR::rasterize(20,
                        operators = "z_p70",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p70.tif")

z_p80 <- lasR::rasterize(20,
                        operators = "z_p80",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p80.tif")

z_p90 <- lasR::rasterize(20,
                        operators = "z_p90",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p90.tif")

z_p95 <- lasR::rasterize(20,
                        operators = "z_p95",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p95.tif")

z_p99 <- lasR::rasterize(20,
                        operators = "z_p99",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_p99.tif")

z_max <- lasR::rasterize(20,
                        operators = "z_max",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_max.tif")

z_mean <- lasR::rasterize(20,
                        operators = "z_mean",
                        filter = c(keep_first(), keep_z_above(1.3)),
                        ofile = "./metrics/lidar/z_mean.tif")

z_sd <- lasR::rasterize(20,
                         operators = "z_sd",
                         filter = c(keep_first(), keep_z_above(1.3)),
                         ofile = "./metrics/lidar/z_sd.tif")

z_cv <- lasR::rasterize(20,
                          operators = "z_cv",
                          filter = c(keep_first(), keep_z_above(1.3)),
                          ofile = "./metrics/lidar/z_cv.tif")

z_above2 <- lasR::rasterize(20,
                            operators = "z_above2",
                            filter = keep_first(),
                            ofile = "./metrics/lidar/z_above2.tif")

z_above5 <- lasR::rasterize(20,
                            operators = "z_above5",
                            filter = keep_first(),
                            ofile = "./metrics/lidar/z_above5.tif")

z_above10 <- lasR::rasterize(20,
                             operators = "z_above10",
                             filter = keep_first(),
                             ofile = "./metrics/lidar/z_above10.tif")

z_above15 <- lasR::rasterize(20,
                             operators = "z_above15",
                             filter = keep_first(),
                             ofile = "./metrics/lidar/z_above15.tif")

z_kurt <- lasR::rasterize(20,
                          operators = "z_kurt",
                          ofile = "./metrics/lidar/z_kurt.tif")

z_skew <- lasR::rasterize(20,
                          operators = "z_skew",
                          ofile = "./metrics/lidar/z_skew.tif")

# Cannot use parallel computing
# # Fractional cover 0.5 m to 2 m
# fractional_cover_05_2 <- lasR::rasterize(20,
#                                          fractional.cover(Z, 0.5, 2),
#                                          filter = "Z < 2",
#                                          "./metrics/lidar/fractional_cover_05_2.tif")
#
# # Fractional cover 2 m to 5 m
# fractional_cover_2_5 <- lasR::rasterize(20,
#                                         fractional.cover(Z, 2, 5),
#                                         filter = "Z < 5",
#                                         "./metrics/lidar/fractional_cover_2_5.tif")
#
# # Fractional cover 5 m to 10 m
# fractional_cover_5_10 <- lasR::rasterize(20,
#                                          fractional.cover(Z, 5, 10),
#                                          filter = "Z < 10",
#                                          "./metrics/lidar/fractional_cover_5_10.tif")
#
# pipeline <- read + normalize + fractional_cover_05_2 + fractional_cover_2_5 + fractional_cover_5_10

# Count of data lower than 0.5
count_0_05 <- lasR::rasterize(20,
                              "count",
                              filter = "Z < 0.5",
                              "./metrics/lidar/count_0_05.tif")

# Count of data lower than 2
count_0_2 <- lasR::rasterize(20,
                             "count",
                             filter = "Z < 2",
                             "./metrics/lidar/count_0_2.tif")

# Count of data lower than 5
count_0_5 <- lasR::rasterize(20,
                             "count",
                             filter = "Z < 5",
                             "./metrics/lidar/count_0_5.tif")

# Count of data lower than 10
count_0_10 <- lasR::rasterize(20,
                              "count",
                              filter = "Z < 10",
                              "./metrics/lidar/count_0_10.tif")





# 🟡 Apply pipe only for necessary metrics ----
# dem
lasR::exec(read + dem,
           on = ctg,
           progress = TRUE,
           ncores = 20)

# veg
lasR::exec(read + normalize + z_p80 + z_p95 + z_cv + z_above2 + count_0_05 + count_0_2 + z_skew,
           on = ctg,
           progress = TRUE,
           ncores = nested(15,30))






# 🟡 Apply pipe ----
# dem
lasR::exec(read + dem,
           on = ctg,
           progress = TRUE,
           ncores = 20)

# z_p
lasR::exec(read + normalize + z_p5 + z_p10 + z_p20 + z_p30 + z_p40 + z_p50 + z_p60 + z_p70 + z_p80 + z_p90 + z_p95 + z_p99,
           on = ctg,
           progress = TRUE,
           ncores = 20)

# z_dist
lasR::exec(read + normalize + z_max + z_mean + z_sd + z_cv,
           on = ctg,
           progress = TRUE,
           ncores = 20)

# z_above
lasR::exec(read + normalize + z_above2 + z_above5 + z_above10 + z_above15,
           on = ctg,
           progress = TRUE,
           ncores = 20)

# z_struct
lasR::exec(read + normalize + z_kurt + z_skew,
           on = ctg,
           progress = TRUE,
           ncores = 20)

# z_count
lasR::exec(read + normalize + count_0_05 + count_0_2 + count_0_5 + count_0_10,
           on = ctg,
           progress = TRUE,
           ncores = 20)





# 🟡 Compute metrics from lidR ----
# Rumple index
plan(multisession, workers = 5L)
set_lidr_threads(0)

# Reduce catalog size based on what have been computed
list.files("D:/00_Ontario_eFRI/OVF/metrics/lidar/rumple_index_temp", pattern = "\\.tif$") %>%
  gsub("_rumple_index.tif", ".laz", .) -> computed_files

ctg %>%
  st_as_sf() %>%
  mutate(filename_short = gsub("D:\\\\00_Ontario_eFRI\\\\OVF\\\\laz\\\\", "", filename)) %>%
  filter(!filename_short %in% computed_files) %>%
  st_buffer(100) %>%
  st_filter(ctg %>% st_as_sf(), .) %>%
  pull(filename) %>%
  readLAScatalog() -> ctg

opt_output_files(ctg) <- "D:/00_Ontario_eFRI/OVF/metrics/lidar/rumple_index_temp/{*}_rumple_index"

pixel_metrics(ctg, ~lidR::rumple_index(x = X,
                                       y = Y,
                                       z = Z),
              res = 20)

list.files("D:/00_Ontario_eFRI/OVF/metrics/lidar/rumple_index_temp", pattern = "\\.tif$", full.names = TRUE) %>%
  vrt("D:/00_Ontario_eFRI/OVF/metrics/lidar/rumple_index_temp/rumple_index.vrt")

rast("D:/00_Ontario_eFRI/OVF/metrics/lidar/rumple_index_temp/rumple_index.vrt") %>%
  writeRaster("D:/00_Ontario_eFRI/OVF/metrics/lidar/rumple_index.tif")





# 🟡 Fill NA with 0 for metrics with "keep_z_above(1.3)" ----
# Reference metric with no filter
rast("./metrics/lidar/z_skew.tif") -> ref

for(m in c("z_cv", "z_p80", "z_p95")){
  print(m)
  rast(paste0("./metrics/lidar/temp/", m, "_raw.tif")) -> metric

  ifel(is.na(metric) & !is.na(ref), 0, metric) -> metric_fill

  metric_fill %>%
    writeRaster(paste0("./metrics/lidar/", m, ".tif"))
}








# 🟡 Metrics from others rasters ----
# Fractional cover 0.5 m to 2 m
rast("./metrics/lidar/temp/count_0_05.tif") -> count_0_05
rast("./metrics/lidar/temp/count_0_2.tif") -> count_0_2
(count_0_2-count_0_05)/count_0_2 -> fractional_cover_05_2

fractional_cover_05_2 %>%
  writeRaster("./metrics/lidar/fractional_cover_05_2.tif")


# Fractional cover 2 m to 5 m
rast("./metrics/lidar/temp/count_0_2.tif") -> count_0_2
rast("./metrics/lidar/temp/count_0_5.tif") -> count_0_5
(count_0_5-count_0_2)/count_0_5 -> fractional_cover_2_5

fractional_cover_2_5 %>%
  writeRaster("./metrics/lidar/fractional_cover_2_5.tif")


# Fractional cover 5 m to 10 m
rast("./metrics/lidar/temp/count_0_5.tif") -> count_0_5
rast("./metrics/lidar/temp/count_0_10.tif") -> count_0_10
(count_0_10-count_0_5)/count_0_10 -> fractional_cover_5_10

fractional_cover_5_10 %>%
  writeRaster("./metrics/lidar/fractional_cover_5_10.tif")


# avg/p99
rast("./metrics/lidar/temp/RMF_20m_T130cm_avg.tif") -> avg
rast("./metrics/lidar/temp/RMF_20m_T130cm_p99.tif") -> p99

avg/p99 -> avg_p99

ifel(!is.nan(avg) & is.nan(avg_p99), 1, avg_p99) -> avg_p99

avg_p99 %>%
  writeRaster("./metrics/lidar/avg_p99.tif")


# Dem from .vrt
# 20 m
rast("E:/dtm_OVF/dtm.vrt") %>%
  terra::aggregate(40, "mean") -> dem

dem %>%
  writeRaster("./metrics/lidar/dem.tif")

# 5 m (OVF)
rast("E:/dtm_OVF/dtm.vrt") %>%
  terra::aggregate(10, "mean") -> dem_5m

# 5 m (RMF)
rast("K:/RMF/RMF_SPL/Romeo_Malette_lidar_drive/DEM/dtm.vrt") %>%
  terra::aggregate(10, "mean") -> dem_5m

dem_5m %>%
  writeRaster("./metrics/lidar/dem_5m.tif")


# DEM at 1 m resolution over drone flight area based on provincial lidar data on metafor-5
st_read("D:/00_Ontario_eFRI/data/drone_flight_areas/drone_flight_areas.shp") %>%
  st_transform(st_crs(ctg)) -> drone_flight_areas

ctg@data %>%
  st_geometry() %>%
  st_union() %>%
  st_filter(drone_flight_areas, .) %>%
  st_buffer(2000) -> drone_flight_areas

map(drone_flight_areas$sector,
    function(x){
      drone_flight_areas %>%
        filter(sector == x) %>%
        st_filter(ctg@data, .) %>%
        pull(filename) %>%
        readLAScatalog() -> ctg_clip

      lasR::exec(lasR::reader() + lasR::dtm(1),
                 on = ctg_clip,
                 progress = TRUE,
                 ncores = 10) -> dem

      dem %>%
        writeRaster(paste0("./metrics/lidar/dem_", x, "_1m.tif"))
    })



# Slope
rast("./metrics/lidar/dem.tif") -> dem
dem %>%
  terrain(v = "slope",
          unit = "radians") %>%
  tan() %>%
  {.*100} -> slope

slope %>%
  writeRaster("./metrics/lidar/slope.tif")


# Aspect
rast("./metrics/lidar/dem.tif") -> dem
dem %>%
  terra::terrain(v = "aspect",
                 unit = "degrees") -> aspect

aspect %>%
  writeRaster("./metrics/lidar/aspect.tif")


# Relative topographic position
wbt_relative_topographic_position(dem = "./metrics/lidar/dem.tif",
                                  output = "./metrics/lidar/relative_topographic_position_60m.tif",
                                  filterx = 3,
                                  filtery = 3)

wbt_relative_topographic_position(dem = "./metrics/lidar/dem.tif",
                                  output = "./metrics/lidar/relative_topographic_position_100m.tif",
                                  filterx = 5,
                                  filtery = 5)

wbt_relative_topographic_position(dem = "./metrics/lidar/dem.tif",
                                  output = "./metrics/lidar/relative_topographic_position_200m.tif",
                                  filterx = 10,
                                  filtery = 10)

wbt_relative_topographic_position(dem = "./metrics/lidar/dem.tif",
                                  output = "./metrics/lidar/relative_topographic_position_400m.tif",
                                  filterx = 20,
                                  filtery = 20)

wbt_relative_topographic_position(dem = "./metrics/lidar/dem.tif",
                                  output = "./metrics/lidar/relative_topographic_position_800m.tif",
                                  filterx = 60,
                                  filtery = 60)


# DEM breached
# 20 m
wbt_breach_depressions(dem = "./metrics/lidar/dem.tif",
                       output = "./metrics/lidar/temp/dem_breached.tif",
                       flat_increment = 0.0001,
                       fill_pits = TRUE)

# 5 m
wbt_breach_depressions(dem = "./metrics/lidar/dem_5m.tif",
                       output = "./metrics/lidar/temp/dem_breached_5m.tif",
                       flat_increment = 0.0001,
                       fill_pits = TRUE)


# SAGA wetness index
# 20 m
RSAGA::rsaga.wetness.index(in.dem = "./metrics/lidar/temp/dem_breached.tif",
                           out.wetness.index = "./metrics/lidar/sagawi.tif",
                           suction = 2000,
                           area.type = "absolute",
                           slope.type = "local",
                           env = env_saga)

# 5 m
RSAGA::rsaga.wetness.index(in.dem = "./metrics/lidar/temp/dem_breached_5m.tif",
                           out.wetness.index = "./metrics/lidar/sagawi_5m.tif",
                           suction = 2000,
                           area.type = "absolute",
                           slope.type = "local",
                           env = env_saga)


# Topographic wetness index
# 20 m
wbt_d8_flow_accumulation(input = "./metrics/lidar/temp/dem_breached.tif",
                         output = "./metrics/lidar/temp/sca.tif",
                         out_type = "specific contributing area")

wbt_slope(dem = "./metrics/lidar/temp/dem_breached.tif",
          output = "./metrics/lidar/temp/slope.tif",
          units = "degrees")

wbt_wetness_index(sca = "./metrics/lidar/temp/sca.tif",
                  slope = "./metrics/lidar/temp/slope.tif",
                  output = "./metrics/lidar/twi.tif")

file.remove("./metrics/lidar/temp/sca.tif",
            "./metrics/lidar/temp/slope.tif")

# 5 m
wbt_d8_flow_accumulation(input = "./metrics/lidar/temp/dem_breached_5m.tif",
                         output = "./metrics/lidar/temp/sca.tif",
                         out_type = "specific contributing area")

wbt_slope(dem = "./metrics/lidar/temp/dem_breached_5m.tif",
          output = "./metrics/lidar/temp/slope.tif",
          units = "degrees")

wbt_wetness_index(sca = "./metrics/lidar/temp/sca.tif",
                  slope = "./metrics/lidar/temp/slope.tif",
                  output = "./metrics/lidar/twi_5m.tif")

file.remove("./metrics/lidar/temp/sca.tif",
            "./metrics/lidar/temp/slope.tif")
