for(p in c(1, seq(5, 95, 5), 99)){
  assign(paste0("lasR_z_p", p),
         lasR::rasterize(c(1, 20),
                         operators = paste0("z_p", p),
                         ofile = paste0("./z_p", p, ".tif")))
}

lasR_z_cv <- lasR::rasterize(c(1, 20),
                             operators = "z_cv",
                             ofile = "./z_cv.tif")

lasR_z_skew <- lasR::rasterize(c(1, 20),
                               operators = "z_skew",
                               ofile = "./z_skew.tif")

lasR_z_above2 <- lasR::rasterize(c(1, 20),
                                 operators = "z_above2",
                                 ofile = "./z_above2.tif")

lasR_z_above5 <- lasR::rasterize(c(1, 20),
                                 operators = "z_above5",
                                 ofile = "./z_above5.tif")
