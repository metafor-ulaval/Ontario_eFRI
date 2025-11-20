# 🟡🟡 Load libraries🟡🟡 ----
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyFiles)
library(leaflet)
library(leaflet.extras)
library(tidyverse)
library(magrittr)
library(rmapshaper)
library(sf)
library(terra)
#library(knitr)
library(exactextractr)
# remotes::install_github("metafor-ulaval/eFRItools")
library(eFRItools)





# 🟡🟡 Functions 🟡🟡 ----
# 🌐 Fonction qui permet de zoomer sur le shapefile 🌐
set.view.auto <- function(map,
                          shp) {
  # Extraction de l'etendue du shapefile
  bbox <- st_bbox(shp)

  # Calculer la largeur et la hauteur
  width <- bbox["xmax"] - bbox["xmin"]
  height <- bbox["ymax"] - bbox["ymin"]

  # Calculer la taille maximale
  size <- max(width, height)

  # Application de la fonction de régression log pour le niveau de zoom
  slope <- -1.436848  # Coefficient de régression
  intercept <- 8.347449  # Intercept

  # Calculer le niveau de zoom avec la fonction
  zoom_fct <- intercept + slope * log(size)

  # Limiter le niveau de zoom au valeur minimales et maximales de leaflet
  zoom_fct <- max(min(zoom_fct, 18), 2)

  setView(map,
          lng = mean(st_coordinates(shp)[, 1]),
          lat = mean(st_coordinates(shp)[, 2]),
          zoom = zoom_fct) -> zoomed_map

  return(zoomed_map)
}





# 🟡🟡 Parameters 🟡🟡 ----
otb_dir <- "D:/00_Ontario_eFRI/logiciels/OTB-9.1.0-Win64/bin"
Sys.setenv(OTB_MAX_RAM_HINT = "65536") # 64 GO
Sys.setenv(OTB_MEMORY_AVAILABLE = "65536") # 64 GO
Sys.setenv(GDAL_CACHEMAX = "65536") # 64 GO






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
                            uiOutput("segmentation_metrics_list"), # Peut-être mettre un min et un max a cocher
                            uiOutput("dendrometric_list"), # Peut-être mettre un min et un max a cocher
                            checkboxGroupInput("masks",
                                               "Choose the masks to peform segmentation :",
                                               choices = c("Roads" = "roads",
                                                           "Waterbodies" = "waterbodies",
                                                           "Watercourses" = "watercourses"),
                                               selected = c("roads", "waterbodies")),
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
                            textInput("epsg", "Choose a crs :", "EPSG:2958"),
                            shinyDirButton("wd", "Choose a working directory", "Please select a folder"),
                            verbatimTextOutput("selected_wd"),
                            textInput("name", "Choose output name :", ""),
                            actionButton("run", "Compute enhanced forest resources inventory polygons")))
      )
    ),

    # 🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵

    # 🔵 Interface principale 🔵
    mainPanel(
      # Utilisation de fluidRow et column pour ajuster la carte
      fluidRow(
        column(12,
               numericInput("latitude", "Latitude", 46.78831),
               numericInput("longitude", "Longitude", -71.31259),
               sliderInput("distance", "Analysis radius (km)",
                           min = 0, max = 200, value = 25),
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

  # 🟠 Reactive working directory 🟠
  selected_wd_reactive <- reactive({
    req(input$wd)
    parseDirPath(roots, input$wd)
  })

  # Directory displayed in text
  output$selected_wd <- renderText({
    selected_wd_reactive()
  })

  # 🟠 Reactive ctg 🟠
  ctg_area_reactive <- reactive({
    req(selected_wd_reactive())
    st_read(paste0(selected_wd_reactive(), "/ctg/ctg.shp")) %>%
      st_buffer(1) %>%
      st_union() %>%
      st_transform(4326) -> ctg_area

    ctg_area
  })

  # 🟠 Reactive analysis radius 🟠
  cercle_reactive <- reactive({
    req(input$longitude, input$latitude, input$distance)
    st_sfc(st_point(c(input$longitude, input$latitude)), crs = 4326) %>%
      st_buffer(input$distance*1000) -> cercle

    cercle
  })

  # 🟢 Map initialisation 🟢
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addProviderTiles(providers$Esri.WorldImagery)
  })

  # 🟠 Update ctg 🟠
  observe({
    req(ctg_area_reactive())
    leafletProxy("map") %>%
      set.view.auto(ctg_area_reactive()) %>%
      addPolygons(data = ctg_area_reactive())
  })

  # 🟠 Observer les coordonnées du clic sur la carte pour centrer l'affichage des parois 🟠
  observeEvent(input$map_click, {
    click <- input$map_click
    updateNumericInput(session, "latitude", value = click$lat)
    updateNumericInput(session, "longitude", value = click$lng)
  })

  # 🟠 Update circle 🟠
  observe({
    req(cercle_reactive())

    leafletProxy("map") %>%
      clearGroup("recherche") %>%
      set.view.auto(cercle_reactive()) %>%
      addCircleMarkers(lng = input$longitude,
                       lat = input$latitude,
                       group = "recherche",
                       color = "turquoise",
                       radius = 0.5,
                       weight = 0.5,
                       fillOpacity = 1) %>%
      addPolygons(data = cercle_reactive(),
                  color = "turquoise",
                  weight = 2,
                  fillOpacity = 0,
                  group = "recherche",
                  options = pathOptions(interactive = FALSE))
  })

  # 🟠 Read metrics informations 🟠
  metrics_infos_reactive <- reactive({
    req(selected_wd_reactive())
    list.files(paste0(selected_wd_reactive(), "/metrics"), full.names = TRUE) %>%
      read_metrics() -> metrics

    metrics
  })

  # 🟣 Metric list 🟣
  output$segmentation_metrics_list <- renderUI({
    req(metrics_infos_reactive())
    checkboxGroupInput(
      inputId = "segmentation_metrics",
      label = "Choose 3 to 8 metrics to peform segmentation :",
      choices = metrics_infos_reactive() %>%
        dplyr::filter(type == "lidar",
                      resolution == 20) %>%
        dplyr::pull(name),
      selected = c("z_p95", "z_cv", "z_above2", "sagawi", "fractional_cover_05_2")
    )
  })

  # 🟣 Dendrometric list 🟣
  output$dendrometric_list <- renderUI({
    req(metrics_infos_reactive())
    checkboxGroupInput(
      inputId = "selected_dendrometrics",
      label = "Choose at least 3 dendrometric characteristics :",
      choices = metrics_infos_reactive() %>%
        dplyr::filter(type == "dendro",
               resolution == 20) %>%
        dplyr::pull(name),
      selected = c("dens", "qmdbh", "Vmerch_ha" )
    )
  })

  # 🟠 Compute enhanced forest resources inventory polygons 🟠
  observeEvent(input$run, {
    # 🟢 Read data 🟢
    showNotification("Read data", type = "message", duration = 15, session = session)

    # Best models for imputation
    read.csv("D:/00_Ontario_eFRI/RMF/results/best_model.csv") -> best_models

    best_models %>%
      pull(knn_vars) %>%
      strsplit(",") %>%
      unlist() %>%
      unique() -> best_models_variables

    # Catalog
    st_read(paste0(selected_wd_reactive(), "/ctg/ctg.shp")) %>%
      st_transform(input$epsg) -> ctg

    # Metrics
    metrics_infos %>%
      filter(name %in% c(best_models_variables,
                         input$segmentation_metrics,
                         input$selected_dendrometrics,
                         "z_p95", "z_above2", "fractional_cover_05_2", "slope", "sagawi")) %T>%
      {pull(.,name) ->> metrics_names} %>% # Extract metrics names in the right order
      pull(path) %>%
      map(rast) %>%
      map(project, input$epsg, method = "bilinear") %>%
      map(resample, .[[1]]) %>%
      rast -> metrics

    # Assign metrics names
    names(metrics) <- metrics_names

    # Masks
    paste0(selected_wd_reactive(), "/shp/", input$masks, ".shp") %>%
      map(vect) %>%
      map(project, input$epsg) -> masks

    # Forest inventory polygons (fri)
    st_read(paste0(selected_wd_reactive(), "/shp/PolygonForest.shp")) %>%
      rowid_to_column("id") %>%
      st_transform(input$epsg) -> fri_polygons

    # Landcover
    rast(paste0(selected_wd_reactive(), "/metrics/other/landcover.tif")) %>%
      terra::project(input$epsg, method = "near") -> landcover

    landcover_codes <- c(`1` = "Clear_Open_Water",
                         `2` = "Turbid_Water",
                         `3` = "Shoreline",
                         `4` = "Mudflats",
                         `5` = "Marsh",
                         `6` = "Swamp",
                         `7` = "Fen",
                         `8` = "Bog",
                         `10` = "Heath",
                         `11` = "Sparse_Treed",
                         `12` = "Treed_Upland",
                         `13` = "Deciduous_Treed",
                         `14` = "Mixed_Treed",
                         `15` = "Coniferous_Treed",
                         `16` = "Plantations_Treed_Cultivated",
                         `17` = "Hedge_Rows",
                         `18` = "Disturbance",
                         `19` = "Open_Cliff_Talus",
                         `20` = "Alvar",
                         `21` = "Sand_Barren_Dune",
                         `22` = "Open_Tallgrass_Prairie",
                         `23` = "Tallgrass_Savannah",
                         `24` = "Tallgrass_Woodland",
                         `25` = "Sand_Gravel_Mine_Tailings_Extraction",
                         `26` = "Bedrock",
                         `27` = "Communit_Infrastructure",
                         `28` = "Agriculture_Undifferentiated_Rural_Land_Use",
                         `157` = "Other",
                         `247` = "Cloud_Shadow")

    # Attach levels (lookup table)
    levels(landcover) <- data.frame(value = as.integer(names(landcover_codes)),
                                    class = unname(landcover_codes))

    # Forest_Age_1985-2020
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_age_2019.tif")) %>%
      terra::project(input$epsg, method = "near") -> forest_age_2019

    # Forest_Fire_1985-2020
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_fire_1985_2020.tif")) %>%
      terra::project(input$epsg, method = "near") -> forest_fire_1985_2020

    # Forest_Harvest_1985-2020
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_harvest_1985_2020.tif")) %>%
      terra::project(input$epsg, method = "near") -> forest_harvest_1985_2020



    # CLIP ALL DATA #
    cercle_reactive() %>%
      st_as_sf() %>%
      st_transform(input$epsg) -> cercle_clip

    ctg %<>%
      ms_clip(cercle_clip)

    metrics %<>%
      crop(cercle_clip) %>%
      mask(cercle_clip)

    masks %<>%
      map(crop, cercle_clip)

    fri_polygons %<>%
      st_filter(cercle_clip) %>%
      ms_clip(cercle_clip)

    landcover %<>%
      crop(cercle_clip) %>%
      mask(cercle_clip)

    forest_age_2019 %<>%
      crop(cercle_clip) %>%
      mask(cercle_clip)

    forest_fire_1985_2020 %<>%
      crop(cercle_clip) %>%
      mask(cercle_clip)

    forest_harvest_1985_2020 %<>%
      crop(cercle_clip) %>%
      mask(cercle_clip)


    # 🟢 Create and set wd 🟢
    segmentation_wd <- paste0(selected_wd_reactive(), "/segmentations/", input$name)
    dir.create(segmentation_wd)

    # 🟢 Save segmentation information 🟢
    cercle_clip %>%
      st_write(paste0(segmentation_wd, "/cercle_clip.shp"))



    # 🟢 Segmentation 🟢
    showNotification("Perform segmentation", type = "message", duration = 15, session = session)

    eFRI_segmentation(metrics = metrics[[input$segmentation_metrics]],
                      masks = masks,
                      thresh = input$grm_thresh,
                      spec = input$grm_spec,
                      spat = input$grm_spat,
                      method = "bs",
                      output_path = segmentation_wd,
                      output_name = "segmentation",
                      otb_dir = otb_dir)

    # 🟢 Imputation 🟢
    showNotification("Peform imputation", type = "message", duration = 15, session = session)

    eFRI_imputation(segmentation = st_read(paste0(segmentation_wd, "/segmentation.shp")),
                    segmentation_id_field = "id",
                    forest_polygon = fri_polygons,
                    metrics = metrics,
                    landcover = landcover,
                    forest_fire = forest_fire_1985_2020,
                    forest_harvest = forest_harvest_1985_2020,
                    ctg = ctg,
                    lidar_year_field = "Fl_Cr_Y",
                    forest_year_field = "YRUPD",
                    forest_composition_field = "SPCOMP",
                    forest_type_field = "POLYTYPE",
                    target_var = best_models$target_var,
                    knn_var = best_models$knn_vars) -> segmentation_imputed

    # 🟢 Build attribute table 🟢
    showNotification("Build attribute table in segmented polygons", type = "message", duration = 15, session = session)

    eFRI_attribute_table(segmentation = segmentation_imputed,
                         metrics = metrics,
                         selected_dendrometrics = input$selected_dendrometrics,
                         landcover = landcover,
                         forest_fire = forest_fire_1985_2020,
                         forest_harvest = forest_harvest_1985_2020,
                         forest_age = forest_age_2019) -> segmentation_data

    segmentation_data %>%
      st_write(dsn = paste0(segmentation_wd, "/segmentation.gpkg"), layer = "segmentation_data")

    showNotification("Done !", type = "message", duration = 15, session = session)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
