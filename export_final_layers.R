setwd("D:/00_Ontario_eFRI/10_livrables")

# RMF
st_read("./segmentation/RMF/model_B/data.gpkg",
        layer = "model_B_data_imputation_smoothed_clipped") %>%
dplyr::select(-InPoly_FID, -SmoPgnFlag) %>%
st_write("./segmentation/RMF/model_B/Automated_stand_delineation_RMF.gpkg",
        layer = "Automated_stand_delineation_RMF")

# OVF
st_read("./segmentation/OVF/model_B/data.gpkg",
        layer = "model_B_data_imputation_smoothed") %>%
  dplyr::select(-InPoly_FID, -SmoPgnFlag) %>%
st_write("./segmentation/OVF/model_B/Automated_stand_delineation_OVF.gpkg",
         layer = "Automated_stand_delineation_OVF")
