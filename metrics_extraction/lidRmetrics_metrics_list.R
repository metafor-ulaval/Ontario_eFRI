metrics_percabove_mean <- function(z) {

  n <- length(z)
  pzabovemean <- lidR:::fast_countover(z, mean(z))/n * 100
  return(pzabovemean)

}

lidRmetrics_z_abovemean <- lasR::rasterize(c(1, 20),
                                           operators = metrics_percabove_mean(Z),
                                           ofile = "./z_abovemean.tif")

lidRmetrics_z_pcum <- lasR::rasterize(c(1, 20),
                                      operators = lidRmetrics::metrics_canopydensity(Z),
                                      ofile = "./z_pcum.tif")

lidRmetrics_pz <- lasR::rasterize(c(1, 20),
                                  operators = lidRmetrics::metrics_interval(Z),
                                  ofile = "./pz.tif")

lidRmetrics_dispersion <- lasR::rasterize(c(1, 20),
                                          operators = lidRmetrics::metrics_dispersion(Z),
                                          ofile = "./dispersion.tif")
