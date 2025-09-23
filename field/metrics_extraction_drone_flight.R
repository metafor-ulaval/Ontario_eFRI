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





# 🟡 Setwd, parameters and functions ----
setwd("F:/Ontario_eFRI/02_donnees_traitees")

plan(multisession, workers = 5L)
set_lidr_threads(0)

list.files() -> sub_sectors





# 🟡 Computing ----
sub_sector <- sub_sectors[[1]]
for(sub_sector in sub_sectors){
  # Read data
  # Catalog
  readLAScatalog(paste0("./", sub_sector ,"/nuages/01_classified")) -> ctg
  readLAScatalog(paste0("./", sub_sector ,"/nuages/02_normalized")) -> ctg_normalized






  # # Compute all lidRmetrics except texture metrics
  # all_lidr_metrics <- pixel_metrics(ctg_normalized, ~lidRmetrics::metrics_set3(x = X,
  #                                                                   y = Y,
  #                                                                   z = Z,
  #                                                                   i = Intensity,
  #                                                                   ReturnNumber = ReturnNumber,
  #                                                                   NumberOfReturns = NumberOfReturns,
  #                                                                   zmin = NA,
  #                                                                   threshold = c(2, 5),
  #                                                                   dz = 1,
  #                                                                   interval_count = 10,
  #                                                                   zintervals = c(0, 0.15, 2, 5, 10, 20, 30),
  #                                                                   pixel_size = 1,
  #                                                                   vox_size = 1,
  #                                                                   KeepReturns = c(1, 2, 3, 4)),
  #                                   res = 20)
  # 
  # # Write individual metrics
  # for(metric in names(all_lidr_metrics)){
  #   cat(paste0(metric, "\n"))
  #   all_lidr_metrics[[metric]] %>%
  #     writeRaster(paste0("./", sub_sector ,"/metrics/", metric, ".tif"))
  # }
  # 
  # 
  # 
  # 
  # 
  # # DEM 1 m and 20 m resolution
  # #read <- reader("-keep_random_fraction 0.1")
  # read <- reader()
  # 
  # dem_20m <- lasR::dtm(20,
  #                      ofile = paste0("./", sub_sector ,"/metrics/dem_20m.tif"))
  # 
  # dem_1m <- lasR::dtm(1,
  #                     ofile = paste0("./", sub_sector ,"/metrics/dem_1m.tif"))
  # 
  # # Apply pipe
  # lasR::exec(read + dem_20m + dem_1m,
  #            on = ctg,
  #            progress = TRUE,
  #            ncores = nested(ncores = ceiling((ncores()-4)/4L), ncores2 = 4L))
  # 
  # 
  # 
  # 
  # 
  # # SAGAWI 1 m resolution
  # wbt_feature_preserving_smoothing(dem = paste0("./", sub_sector ,"/metrics/dem_1m.tif"),
  #                                  output = paste0("./", sub_sector ,"/metrics/dem_smoothed_1m.tif"),
  #                                  filter = 15,
  #                                  num_iter = 3,
  #                                  norm_diff = 30)
  # 
  # wbt_breach_depressions(dem = paste0("./", sub_sector ,"/metrics/dem_smoothed_1m.tif"),
  #                        output = paste0("./", sub_sector ,"/metrics/dem_breached_1m.tif"),
  #                        flat_increment = 0.0001,
  #                        fill_pits = TRUE)
  # 
  # RSAGA::rsaga.wetness.index(in.dem = paste0("./", sub_sector ,"/metrics/dem_breached_1m.tif"),
  #                            out.wetness.index = paste0("./", sub_sector ,"/metrics/sagawi_1m.tif"),
  #                            suction = 10,
  #                            area.type = "absolute",
  #                            slope.type = "local",
  #                            env = env_saga)
  # 
  # 
  # 
  # 

  # Buffered Area Based Approach
    # Create pipe step
    read <- reader(filter = paste0("-keep_random_fraction", 0.5))


    z_p95 <- lasR::rasterize(c(2, 20),
                             operators = "z_p95",
                             filter = c(keep_first(), keep_z_above(1.3)),
                             ofile = paste0("./", sub_sector ,"/metrics/z_p95_MWABA_frac_", frac*100, "_buffer_", w,"m.tif"))

    z_above2 <- lasR::rasterize(c(2, 20),
                                operators = "z_above2",
                                filter = keep_first(),
                                ofile = paste0("./", sub_sector ,"/metrics/z_above2_MWABA_frac_", frac*100, "_buffer_", w,"m.tif"))

    z_cv <- lasR::rasterize(c(2, 20),
                            operators = "z_cv",
                            filter = c(keep_first(), keep_z_above(1.3)),
                            ofile = paste0("./", sub_sector ,"/metrics/z_cv_MWABA_frac_", frac*100, "_buffer_", w,"m.tif"))
    
    start_time <- Sys.time()

    # Apply pipe
    lasR::exec(read + z_p95 + z_above2 + z_cv,
               on = ctg_normalized,
               progress = TRUE,
               ncores = nested(ncores = ceiling((ncores()-4)/8L), ncores2 = 8L))
    
    end_time <- Sys.time()
    
    end_time - start_time

}