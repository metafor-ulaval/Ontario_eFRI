library(tidyverse)
library(sf)

setwd("F:/Ontario_eFRI/02_donnees_traitees")

# 🟡 Computing ----
map_dfr(1:10,~{
  lidR::readLAScatalog(paste0("./S", .x ,"/nuages/02_norm_decimated"))@data %>%
    st_geometry() %>%
    st_buffer(10) %>%
    st_union() %>%
    st_as_sf() %>%
    st_transform(4326) %>%
    mutate(sector = paste0("S", .x))
}) -> data

data %>%
  st_write("drone_flight_areas.shp")
