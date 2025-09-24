# 🟡🟡 Load libraries🟡🟡 ----
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyFiles)
library(tidyverse)
library(magrittr)
library(sf)
library(terra)
library(smoothr)
library(pfif)
library(knitr)
library(exactextractr)
library(leaflet)
library(leaflet.extras)
library(basemapR)





# 🟡🟡 Imputation functions 🟡🟡 ----
rstudioapi::getSourceEditorContext()$path %>% 
  dirname() %>% 
  paste0("/functions") %>% 
  list.files(full.names = T) %>% 
  map(source)






# 🟡🟡 Set wd and paramters 🟡🟡 ----
"D:/00_Ontario_eFRI/RMF" -> wd

epsg <- "EPSG:2958" # RMF
otb_dir <- "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin"
Sys.setenv(OTB_MAX_RAM_HINT = "65536") # 64 GO
Sys.setenv(OTB_MEMORY_AVAILABLE = "65536") # 64 GO
Sys.setenv(GDAL_CACHEMAX = "65536") # 64 GO

best_models <- read.csv(paste0(wd, "/results/best_model.csv"))

best_models %>% 
  pull(knn_vars) %>% 
  strsplit(",") %>% 
  unlist() %>% 
  unique() -> best_models_variables






# 🟡🟡 Read metrics informations 🟡🟡 ----
lidar_metrics_desc <- c("dem" = "digital elevation model",
                        "dem_5m" = "digital elevation model",
                        "slope" = "slope in percent",
                        "aspect" = "slope orientation in degrees",
                        "z_above2" = "canopy cover of first returns above 2 m",
                        "z_above5" = "canopy cover of first returns above 5 m",                        
                        "z_above10" = "canopy cover of first returns above 10 m", 
                        "z_above15" = "canopy cover of first returns above 15 m",
                        "z_mean" = "average height of returns above 1.3 m",
                        "z_cv" = "coefficient of variation of returns above 1.3 m",
                        "z_kurt" = "kurtosis  of returns above 1.3 m",
                        "z_max" = "maximum height of returns above 1.3 ms",
                        "z_p5" = "height of 5th percentiles of returns above 1.3 m",
                        "z_p10" = "height of 10th percentiles of returns above 1.3 m",
                        "z_p20" = "height of 20th percentiles of returns above 1.3 m",
                        "z_p30" = "height of 30th percentiles of returns above 1.3 m",
                        "z_p40" = "height of 40th percentiles of returns above 1.3 m",
                        "z_p50" = "height of 50th percentiles of returns above 1.3 m",
                        "z_p60" = "height of 60th percentiles of returns above 1.3 m",
                        "z_p70" = "height of 70th percentiles of returns above 1.3 m",
                        "z_p80" = "height of 80th percentiles of returns above 1.3 m",
                        "z_p90" = "height of 90th percentiles of returns above 1.3 m",
                        "z_p95" = "height of 95th percentiles of returns above 1.3 m",
                        "z_p99" = "height of 99th percentiles of returns above 1.3 m",
                        "rumple_index" = "rumple index of raw lidar data",
                        "z_qav" = "average square height of returns above 1.3 m",
                        "z_skew" = "skewness of returns above 1.3 m",
                        "z_sd" = "standard deviation of returns above 1.3 m",
                        "avg_p99" = "avg/p99",
                        "fractional_cover_05_2" = "number of returns between 0.5 m and 2 m divided by all returns below 2 m",
                        "fractional_cover_2_5" = "number of returns between 2 m and 5 m divided by all returns below 5 m",
                        "fractional_cover_5_10" = "number of returns between 5 m and 10 m divided by all returns below 10 m",
                        "relative_topographic_position_60m" = "relative topographic position with 3 pixels in x and y (60 m) from whiteboxtools",
                        "relative_topographic_position_100m" = "relative topographic position with 5 pixels in x and y (100 m) from whiteboxtools",
                        "relative_topographic_position_200m" = "relative topographic position with 10 pixels in x and y (200 m) from whiteboxtools",
                        "relative_topographic_position_400m" = "relative topographic position with 20 pixels in x and y (400 m) from whiteboxtools",
                        "relative_topographic_position_800m" = "relative topographic position with 40 pixels in x and y (800 m) from whiteboxtools",
                        "twi" = "topographic wetness index",
                        "twi_5m" = "topographic wetness index",
                        "sagawi" = "saga wetness index",
                        "sagawi_5m" = "saga wetness index")

tibble(path = list.files(paste0(wd, "/metrics/lidar"), pattern = "\\.tif$", full.names = TRUE)) %>% 
  mutate(name = list.files(paste0(wd, "/metrics/lidar"), pattern = "\\.tif$") %>% gsub(".tif", "", .),
         desc = lidar_metrics_desc[name],
         res = path %>% 
           map(rast) %>% 
           map(res) %>% 
           map(1) %>% 
           unlist(),
         type = "lidar") -> lidar_metrics_infos


# sentinel2
sentinel2_metrics_desc <- c("B2" = "Blue - 490 nm",
                            "B3" = "Green - 560 nm",
                            "B4" = "Red - 665 nm",
                            "B5" = "Vegetation Red Edge 1 (RE1) - 705 nm",
                            "B6" = "Vegetation Red Edge 2 (RE2) - 740 nm",
                            "B7sr" = "Surface Reflectance Vegetation Red Edge 3 - 783 nm",
                            "B8" = "Near Infrared (NIR) - 842 nm",
                            "B8sr" = "Surface Reflectance Near Infrared (NIRsr) - 842 nm",
                            "B8Asr" = "Surface Reflectance Narrow (Nsr) - 865 nm",
                            "B11" = "Shortwave Infrared 1 (SWIR1) - 1610 nm",
                            "B12" = "Shortwave Infrared 2 (SWIR2) - 2200 nm",
                            "NDVI" = "Normalized Difference Vegetation Index ((NIR - Red) / (NIR + Red))",
                            "GNDVI" = "Green Normalized Difference Vegetation Index ((NIR - Green) / (NIR + Green))",
                            "NDVIre" = "Normalized Difference Vegetation Index using Red Edge ((RE1 - Red) / (RE1 + Red))",
                            "RVI" = "Ratio Vegetation Index (NIR / Red)",
                            "DVI" = "Difference Vegetation Index (NIR - Red)",
                            "SAVI" = "Soil Adjusted Vegetation Index ((NIR - Red) / (NIR + Red + 0.5))*(1+0.5)",  
                            "EVI" = "Enhanced Vegetation Index 2.5*((NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1))",
                            "GCI" = "Green Chlorophyll Index ((NIR / Green) - 1)",
                            "NDMI" = "Normalized Difference Moisture Index ((NIR - SWIR1) / (NIR + SWIR1))",  
                            "NDWI" = "Normalized Difference Water Index ((NIR - SWIR2) / (NIR + SWIR2))")

tibble(path = list.files(paste0(wd, "/metrics/sentinel2"), pattern = "\\.tif$", full.names = TRUE)) %>% 
  mutate(name = list.files(paste0(wd, "/metrics/sentinel2"), pattern = "\\.tif$") %>% gsub(".tif", "", .),
         desc = sentinel2_metrics_desc[name],
         res = path %>% 
           map(rast) %>% 
           map(res) %>% 
           map(1) %>% 
           unlist(),
         type = "sentinel_2") -> sentinel2_metrics_infos

# dendro
dendro_metrics_desc <- c("AGB_ha" = "above ground biomass normalized per hectare (t C/ha)",
                         "AGB_ha_FOR" = "above ground biomass normalized per hectare (t C/ha) / mask for RMF",
                         "ba_ha" = "basal area / tree cross sectional area (approximated as a circle) at breast height (1.3 m) (m²/ha)",
                         "ba_ha_FOR" = "basal area / tree cross sectional area (approximated as a circle) at breast height (1.3 m) (m²/ha) / mask for RMF",
                         "ba_ha_min20" = "basal area / tree cross sectional area (approximated as a circle) at breast height (1.3 m) (m²/ha) / min 20 ??????", # À élucider le min 20
                         "dens" = "stem density with DBH > 7.1 cm (stems/ha)",
                         "dens_FOR" = "tem density with DBH > 7.1 cm (stems/ha) / mask for RMF",
                         "lor" = "lorey’s height /	average tree height weighted by basal area (m)",
                         "lor_FOR" = "lorey’s height /	average tree height weighted by basal area (m) / mask for RMF",
                         "qmdbh" = "quadratic mean of diameter at breast height (cm)",
                         "qmdbh_FOR" = "quadratic mean of diameter at breast height (cm) / mask for RMF",
                         "top_height" = "maximum height (m)",
                         "top_height_FOR" = "maximum height (m) / mask for RMF",
                         "V_ha" = "total whole stem volume per hectare	(m³/ha)",
                         "V_ha_FOR" = "total whole stem volume per hectare	(m³/ha) / mask for RMF",
                         "Vmerch_ha" = "total merchantable volume normalized per hectare where stump height is set to 0.2 m and minimum to diameter to 10 cm (m³/ha)",
                         "Vmerch_ha_FOR" = "total merchantable volume normalized per hectare where stump height is set to 0.2 m and minimum to diameter to 10 cm (m³/ha) / mask for RMF")

tibble(path = list.files(paste0(wd, "/metrics/dendro"), pattern = "\\.tif$", full.names = TRUE)) %>% 
  mutate(name = list.files(paste0(wd, "/metrics/dendro"), pattern = "\\.tif$") %>% gsub(".tif", "", .), # CA NE MARCHERA PAS POUR OVF, RENOMMER PLUS SIMPLEMENT
         desc = dendro_metrics_desc[name],
         res = path %>% 
           map(rast) %>% 
           map(res) %>% 
           map(1) %>% 
           unlist(),
         type = "dendrometric") %>% 
  filter(!grepl("_FOR$", name)) -> dendro_metrics_infos

bind_rows(lidar_metrics_infos,
          sentinel2_metrics_infos,
          dendro_metrics_infos) -> metrics_infos

rm(dendro_metrics_desc,
   lidar_metrics_desc,
   sentinel2_metrics_desc,
   lidar_metrics_infos,
   sentinel2_metrics_infos,
   dendro_metrics_infos)





# 🟡🟡 Shiny App 🟡🟡 ----
ui <- fluidPage(
  useShinyjs(),  # pour gérer l'interactivité
  titlePanel("Enhanced Forest Resources Inventory"),
  sidebarLayout(
    
    # 🔵 Barre latérale 🔵
    sidebarPanel(
      bsCollapse(
        open = "Create new segmentation",
        bsCollapsePanel("Create new segmentation",
                        div(class = "sidebar-panel",
                            # Parameters
                            uiOutput("metric_list"), # Peut-être mettre un min et un max a cocher
                            uiOutput("dendrometric_list"), # Peut-être mettre un min et un max a cocher                      
                            checkboxGroupInput("masks", 
                                               "Choose the masks to peform segmentation :", 
                                               choices = c("Roads" = "roads",
                                                           "Waterbodies" = "waterbodies",
                                                           "Watercourses" = "watercourses"), 
                                               selected = c("selected_metric")),
                            strong("Choose generic region merging parameters :"),
                            numericInput("grm_thresh", "Threshold", 47),
                            strong("Maximum heterogeneity allowed when merging, controlling how fine or coarse the segmentation is"),
                            "Small threshold → fine segmentation (many small regions).",
                            "Large threshold → coarse segmentation (fewer, larger regions).",
                            numericInput("grm_spec", "Weight of pectral homogeneity", 0.6),
                            strong("Weight given to spectral similarity (pixel values, band means) when deciding whether regions should merge."),
                            "If high → segmentation will mostly respect spectral values.",
                            "If low → spectral similarity matters little, other criteria (shape) dominate.",
                            numericInput("grm_spat", "Weight of patial homogeneity", 0.6),
                            strong("Weight given to shape similarity (compactness and smoothness) when deciding whether regions should merge."),
                            "If high → the algorithm favors compact, smooth regions even if spectral similarity is weaker.",
                            "If low → region boundaries will mostly follow spectral homogeneity.",
                            # Action button to show plots and compute statistics
                            shinyDirButton("wd", "Choose a working directory", "Please select a folder"),
                            verbatimTextOutput("selected_wd"),
                            actionButton("run", "Compute enhanced forest resources inventory polygons"))),
        bsCollapsePanel("Read segmentation",
                        div(class = "sidebar-panel",
                            # Parameters
                            # Make a map to read data with symbology ?
                            ))       
      )
    ),
    
    # 🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵
    
    # 🔵 Interface principale 🔵
    mainPanel(
      # Utilisation de fluidRow et column pour ajuster la carte
      fluidRow(
        column(12,
               leafletOutput("map", height = "600px")  # Carte de taille fixe
        )
      )
    )
    
    # 🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵
    
  )
)





# 🟡🟡 Server 🟡🟡 ----
server <- function(input, output, session) {
  
  # 🟢 = Steps
  # 🟠 = Reactive components
  # 🟣 = Action button components
  
  # 🟠 Give access to the whole filesystem 🟠
  letters <- c("D", "E", "F", "G", "H", "I", "J", "K")
  all_drives <- paste0(letters, ":/")
  roots <- setNames(all_drives[dir.exists(all_drives)],
                    all_drives[dir.exists(all_drives)])
  
  shinyDirChoose(input, "wd", roots = roots)
  
  # Reactive working directory
  selected_wd_reactive <- reactive({
    req(input$wd)
    parseDirPath(roots, input$wd)
  })
  
  # Directory displayed in text
  output$selected_wd <- renderText({
    selected_wd_reactive()
  })
  
  # 🟣 Metric list 🟣
  output$metric_list <- renderUI({
    checkboxGroupInput(
      inputId = "segmentation_metrics",
      label = "Choose 3 to 8 metrics to peform segmentation :",
      choices = metrics_infos %>% 
        filter(type == "lidar",
               res == 20) %>% 
        pull(name),
      selected = c("z_p95", "z_cv", "z_above2", "sagawi", "fractional_cover_05_2")
    )
  })
  
  # 🟣 Dendrometric list 🟣
  output$dendrometric_list <- renderUI({
    checkboxGroupInput(
      inputId = "selected_dendrometrics",
      label = "Choose at least 3 dendrometric characteristics :",
      choices = metrics_infos %>% 
        filter(type == "dendrometric",
               res == 20) %>% 
        pull(name),
      selected = c("dens", "qmdbh", "Vmerch_ha" )
    )
  })
  
  # 🟠 Compute enhanced forest resources inventory polygons 🟠
  observeEvent(input$run, {
    # 🟢 Read data 🟢
    showNotification("Read data", type = "message", duration = 15, session = session)
    # Segmentation data
    # Metrics
    metrics_infos %>%
      filter(name %in% c(best_models_variables, 
                         input$segmentation_metrics, 
                         input$selected_dendrometrics,
                         "z_p95", "z_above2", "fractional_cover_05_2", "slope", "sagawi")) %T>% 
      {pull(.,name) ->> metrics_names} %>% # Extract metrics names in the right order
      pull(path) %>% 
      map(rast) %>% 
      map(project, epsg, method = "bilinear") %>% 
      map(resample, .[[1]]) %>% 
      rast -> metrics
    
    # Assign metrics names
    names(metrics) <- metrics_names
    
    # Masks
    list(vect(paste0(wd, "/shp/roads.shp")),
         vect(paste0(wd, "/shp/waterbodies.shp"))) %>% 
      map(project, epsg) -> masks
    
    # Imputation data
    # Catalog
    st_read(paste0(wd, "/ctg/ctg.shp")) %>% 
      st_transform(epsg) -> ctg
    
    # Subset metrics in a 2 km radius circle in the centroid of area ################### CLIP
    ctg %>% ################### CLIP
      st_bbox() %>% ################### CLIP
      st_as_sfc() %>% ################### CLIP
      st_centroid() %>% ################### CLIP
      st_buffer(1000) %>% ################### CLIP
      st_as_sf() -> subset_circle ################### CLIP
    
    # Forest inventory polygons (fri)
    st_read(paste0(wd, "/shp/PolygonForest.shp")) %>% 
      rowid_to_column("id") %>% 
      st_transform(epsg) -> fri_polygons
    
    # Landcover
    rast(paste0(wd, "/metrics/other/landcover.tif")) -> landcover
    
    # Forest_Age_1985-2020
    rast(paste0(wd, "/metrics/other/forest_age_2019.tif")) -> forest_age_2019
    
    # Forest_Fire_1985-2020
    rast(paste0(wd, "/metrics/other/forest_fire_1985_2020.tif")) -> forest_fire_1985_2020
    
    forest_fire_1985_2020 %>% 
      as.polygons() %>% 
      st_as_sf() %>% 
      dplyr::filter(.[[1]] != 0) %>% 
      st_cast("MULTIPOLYGON") %>% 
      st_cast("POLYGON") %>% 
      rename(year_fire = 1) %>% 
      st_transform(epsg) -> forest_fire_1985_2020_poly
    
    # Forest_Harvest_1985-2020
    rast(paste0(wd, "/metrics/other/forest_harvest_1985_2020.tif")) -> forest_harvest_1985_2020
    
    forest_harvest_1985_2020 %>% 
      as.polygons() %>% 
      st_as_sf() %>% 
      dplyr::filter(.[[1]] != 0) %>% 
      st_cast("MULTIPOLYGON") %>% 
      st_cast("POLYGON") %>% 
      rename(year_harvest = 1) %>% 
      st_transform(epsg) -> forest_harvest_1985_2020_poly
    
    # 🟢 Segmentation 🟢
    showNotification("Perform segmentation", type = "message", duration = 15, session = session)
    ################### CLIP
    metrics %<>% ################### CLIP
      crop(subset_circle) %>% ################### CLIP
      mask(subset_circle) ################### CLIP
    
    # Segmentation
    metrics[[input$segmentation_metrics]] %>%
      pre_processing(masks = masks, smooth = 0) %>% # Normalizes the values to range in ⁠[0,255] with 1 and 99 percentile⁠
      generic_region_merging(thresh = input$grm_thresh,
                             spec = input$grm_spec,
                             spat = input$grm_spat, 
                             simplify = FALSE,
                             ofile = paste0(selected_wd_reactive(), "/grm_raster.tif"), # Save the segmentation in raster
                             otb_dir = otb_dir) -> grm_polygons
    
    # Multipart to single part polygons
    grm_polygons %<>%
      st_cast("MULTIPOLYGON") %>% 
      st_cast("POLYGON") %>% 
      rename("id_segmentation" = 1) %>% 
      rowid_to_column("id") %>% # Compute id
      st_transform(epsg)
    
    # Write polygons
    grm_polygons %>%       
      st_write(paste0(selected_wd_reactive(), "/grm_polygons.gpkg"), 
               layer = "grm_polygons",
               quiet = T)
    
    # 🟢 Smoothing 🟢
    showNotification("Smooth segmentation", type = "message", duration = 15, session = session)
    grm_polygons %>% 
      smoothr::smooth(method = "ksmooth", smoothness = 1) %>%
      st_make_valid() %>% 
      st_write(paste0(selected_wd_reactive(), "/grm_polygons.gpkg"), 
               layer = "grm_polygons_smoothed",
               quiet = T)
    
    # 🟢 Compute values in grm polygons 🟢
    showNotification("Compute values in segmented polygons", type = "message", duration = 15, session = session)
    grm_polygons %<>% 
      mutate(SOURCE = "eFRI") %>% # Source of data
      mutate(YRSOURCE = 2025) %>% # Year of the source
      mutate(area = st_area(.) %>% as.numeric()) %>% 
      mutate(perimeter = st_perimeter(.) %>% as.numeric()) %>% 
      mutate.metrics(metrics, fun = "median") %>% 
      mutate.centroid() %>% 
      mutate.landcover(landcover) %>% 
      mutate.simplify.landcover() %>% # Extract the landcover type and its proportion that cover the largest value of the polygon
      mutate.prop.forested() %>% 
      mutate.proportion(forest_fire_1985_2020, "fire") %>% 
      mutate.simplify.proportion("FIRE", 20) %>% # Extract the fire with the highest proportion : age (YRFIRE) and proportion (FIREPROP) (that cover at least 20 %)  
      mutate.proportion(forest_harvest_1985_2020, "harvest") %>% 
      mutate.simplify.proportion("HARVEST", 20) %>% # Extract the harvest with the highest proportion : age (YRHARVEST) and proportion (HARVESTPROP) (that cover at least 20 %)  
      mutate.age(forest_age_2019, 2019, 2025, "median")
    
    # Extract all columns name that have a proportion of "HARVEST" or "FIRE"
    grm_polygons %>% 
      st_drop_geometry() %>% 
      select(contains("HARVEST_") | contains("FIRE_")) %>% 
      names -> perturbation_column
    
    # Mutate final values
    grm_polygons %<>%
      mutate.origin.age(perturbation_column, 80) %>% # Extract the origin age of the stand (age of the last perturbation that cover more than 80 %)
      mutate(YRORG = ifelse(is.na(YRORG), (2025-age_median_2025), YRORG)) %>% # Assing the origin age if there was no perturbation
      mutate(AGE = 2025-YRORG) %>% # Recalculate the age based on the last perturbation that cover more than 80 %
      mutate.perturbation(perturbation_column, 20) %>% # Extrat the last perburtation with : age (YRDEP), type (DEPTYPE) and proportion (DEPPROP) (that cover at least 20 %)
      mutate(fractional_cover_05_2 = fractional_cover_05_2*100)
    
    # 🟢 Filter FRI polygon for imputation 🟢
    showNotification("Filter FRI polygon for imputation", type = "message", duration = 15, session = session)
    ################### CLIP
    fri_polygons %<>% ################### CLIP
      st_filter(subset_circle) ################### CLIP
    
    # Find ids that had an harvest between inventory and lidar flight so that the inventory is not good
    fri_polygons %>%
      dplyr::select(id, 
                    year = YRUPD) %>% 
      st_intersection(ctg %>% dplyr::select(year_lidar = Fl_Cr_Y)) %>% 
      st_intersection(forest_harvest_1985_2020_poly %>% dplyr::select(year_harvest)) %>% 
      mutate(drop = ifelse(year_harvest > year & year_harvest < year_lidar, 1, 0)) %>% 
      filter(drop == 1) %>% 
      pull(id) %>% 
      unique -> drop_id_harvest
    
    # Find ids that had a fire between inventory and lidar flight so that the inventory is not good
    fri_polygons %>%
      dplyr::select(id, 
                    year = YRUPD) %>% 
      st_intersection(ctg %>% dplyr::select(year_lidar = Fl_Cr_Y)) %>% 
      st_intersection(forest_fire_1985_2020_poly %>% dplyr::select(year_fire)) %>% 
      mutate(drop = ifelse(year_fire > year & year_fire < year_lidar, 1, 0)) %>% 
      filter(drop == 1) %>% 
      pull(id) %>% 
      unique -> drop_id_fire
    
    # Remove polygons based on id
    fri_polygons %<>% 
      filter(!id %in% c(drop_id_harvest, drop_id_fire))
    
    # Extract metrics and variables for filtering
    fri_polygons %<>% 
      mutate.metrics(metrics, fun = "median") %>%
      mutate.centroid() %>% 
      mutate.landcover(landcover) %>%  
      mutate.prop.forested()
    
    # Keep forested poygons based on 3 data source
    fri_polygons %<>% 
      filter(POLYTYPE == "FOR") %>% # keep polygon that are forested based on POLYTYPE
      filter(PROPFORESTED >= 50) %>% # keep polygon that have at least 50 % of forested area
      filter(z_p95 >= 5, z_above2 >= 50) # keep polygon with only high vegetation according to lidar data
    
    # 🟢 Imputation 🟢 
    showNotification("Peform imputation", type = "message", duration = 15, session = session)
    # Compute imputation variables
    fri_polygons %<>% 
      mutate.species.prop("SPCOMP") %>%
      mutate.species.order("SPCOMP") %>%
      mutate.forest.type %>% 
      mutate.functionals.groups
    
    for(target_variable_selected in c("SP_NO_1", "SP_NO_2", "FUNCTIONAL_GROUP_3", "FUNCTIONAL_GROUP_5")){
      
      cat(paste0("Imputation of variable : ", target_variable_selected, "\n"))
      
      # Extract knn variables from the best model
      best_models %>% 
        filter(target_var == target_variable_selected) %>% 
        pull(knn_vars) %>%
        strsplit(",") %>% 
        unlist() -> knn_variables
      
      # Select only target columns of the fri polygons with values (drop_na)
      fri_polygons %>%
        dplyr::select(any_of(c("id", knn_variables, target_variable_selected))) %>% # The unique id, all variables used for knn imputation and the current target variable selected
        drop_na(everything()) -> fri_polygons_selected
      
      # Select only target columns of the grm polygons with values (drop_na)
      grm_polygons %>%
        dplyr::select(any_of(c("id", knn_variables))) %>% # The unique id and all variables used for knn imputation
        drop_na(everything()) -> grm_polygons_selected
      
      # Perform knn imputation
      knn.inputation(reference_polygons = fri_polygons_selected, # With all fri polygons with data (non-na)
                     target_polygons = grm_polygons_selected, # For all grm polygons with data (non-na)
                     knn_variables = knn_variables, # knn variables from the best model
                     target_variables = target_variable_selected, # For the current target variable selected
                     k = 5) %>%
        bind_cols(grm_polygons_selected, .) -> grm_polygons_selected
      
      # A left join is necessary to keep all polygons from original grm data (because some have na)
      grm_polygons_selected %>% 
        st_drop_geometry() %>% 
        dplyr::select(any_of(c("id", 
                               target_variable_selected))) %>%
        left_join(grm_polygons, ., by = "id") -> grm_polygons
    }
    
    # 🟢 Final selection 🟢 
    showNotification("Final selection", type = "message", duration = 15, session = session)
    # Order values ----
    grm_polygons %>%
      dplyr::select(SOURCE,
                    YRSOURCE,
                    AREA = area,
                    PERIMETER = perimeter,
                    PROPFORESTED,
                    POLYTYPE,
                    POLYTYPEPROP,
                    YRORG,
                    AGE,
                    YRDEP,
                    DEPTYPE,
                    DEPPROP,                  
                    YRHARVEST,
                    HARVESTPROP,
                    YRFIRE,
                    FIREPROP,
                    input$selected_dendrometrics,
                    HEIGHT = z_p95,
                    CANOPY_COVER = z_above2,
                    DENSITY = fractional_cover_05_2,
                    SLOPE = slope,
                    MOISTURE = sagawi,
                    LEADSP = SP_NO_1,
                    SECSP = SP_NO_2,
                    FUNCTIONAL_GROUP_3, 
                    FUNCTIONAL_GROUP_5) %>% 
      rename_with(toupper, input$selected_dendrometrics) %>%
      st_write(paste0(selected_wd_reactive(), "/grm_polygons.gpkg"), 
               layer = "grm_polygons_data",
               quiet = T)
    showNotification("Done !", type = "message", duration = 15, session = session)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)