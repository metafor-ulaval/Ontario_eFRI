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
                            uiOutput("summary_metrics_list"), # Peut-être mettre un min et un max a cocher
                            uiOutput("masks_list"),
                            strong("Choose generic region merging parameters :"),
                            numericInput("grm_thresh", "Threshold", 47),
                            # strong("Maximum heterogeneity allowed when merging, controlling how fine or coarse the segmentation is"),
                            # "Small threshold → fine segmentation (many small regions).",
                            # "Large threshold → coarse segmentation (fewer, larger regions).",
                            numericInput("grm_spec", "Weight of spectral homogeneity", 0.6),
                            # strong("Weight given to spectral similarity (pixel values, band means) when deciding whether regions should merge."),
                            # "If high → segmentation will mostly respect spectral values.",
                            # "If low → spectral similarity matters little, other criteria (shape) dominate.",
                            numericInput("grm_spat", "Weight of spatial homogeneity", 0.6),
                            # strong("Weight given to shape similarity (compactness and smoothness) when deciding whether regions should merge."),
                            # "If high → the algorithm favors compact, smooth regions even if spectral similarity is weaker.",
                            # "If low → region boundaries will mostly follow spectral homogeneity.",
                            # Action button to show plots and compute statistics
                            shinyDirButton("wd", "Choose a working directory", "Please select a folder"),
                            verbatimTextOutput("selected_wd"),
                            textInput("name", "Choose output name :", ""),
                            actionButton("run", "Compute enhanced forest resources inventory polygons"),
                            actionButton("add_segmentation", "Add enhanced forest resources inventory polygons to map")))
      )
    ),

    # 🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵

    # 🔵 Interface principale 🔵
    mainPanel(
      # Utilisation de fluidRow et column pour ajuster la carte
      fluidRow(
        column(12,
               numericInput("latitude", "Latitude", value = 46.78831),
               numericInput("longitude", "Longitude", value = -71.31259),
               numericInput("distance", "Analysis radius (km)", value = 25),
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
  metrics_area_reactive <- reactive({
    # req(metrics_infos_reactive())
    #
    # metrics_infos_reactive() %>%
    #   slice(1) %>%
    #   pull(path) %>%
    #   rast() -> raster_temp
    #
    # raster_temp[!is.na(raster_temp)] <- 1
    #
    # raster_temp -> raster_temp
    #
    # raster_temp %>%
    #   as.polygons() %>%
    #   st_as_sf() %>%
    #   st_transform(4326) -> metrics_area

    req(selected_wd_reactive())
    st_read(paste0(selected_wd_reactive(), "/ctg/ctg.shp")) %>%
      st_buffer(1) %>%
      st_union() %>%
      st_transform(4326) -> metrics_area

    metrics_area
  })

  # 🟠 Reactive analysis radius 🟠
  extraction_circle_reactive <- reactive({
    req(input$longitude, input$latitude, input$distance)
    st_sfc(st_point(c(input$longitude, input$latitude)), crs = 4326) %>%
      st_buffer(input$distance*1000) -> extraction_circle

    extraction_circle
  })

  # 🟢 Map initialisation 🟢
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addProviderTiles(providers$Esri.WorldImagery)
  })

  # 🟠 Update ctg and extraction_circle position 🟠
  observe({
    req(metrics_area_reactive())
    leafletProxy("map") %>%
      set.view.auto(metrics_area_reactive()) %>%
      addPolygons(data = metrics_area_reactive(),
                  color = "white",
                  weight = 2,
                  fillOpacity = 0)

    metrics_area_reactive() %>%
      st_transform(4326) %>%
      st_centroid() %>%
      st_coordinates() -> ctg_coordinates
    updateNumericInput(session, "longitude", value = ctg_coordinates[1])
    updateNumericInput(session, "latitude", value = ctg_coordinates[2])
  })

  # 🟠 Observer les coordonnées du clic sur la carte pour centrer l'affichage des parois 🟠
  observeEvent(input$map_click, {
    click <- input$map_click
    updateNumericInput(session, "longitude", value = click$lng)
    updateNumericInput(session, "latitude", value = click$lat)
  })

  # 🟠 Update circle 🟠
  observe({
    req(extraction_circle_reactive())

    leafletProxy("map") %>%
      clearGroup("recherche") %>%
      set.view.auto(extraction_circle_reactive()) %>%
      addCircleMarkers(lng = input$longitude,
                       lat = input$latitude,
                       group = "recherche",
                       color = "turquoise",
                       radius = 0.5,
                       weight = 0.5,
                       fillOpacity = 1) %>%
      addPolygons(data = extraction_circle_reactive(),
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
    selectizeInput(
      inputId = "segmentation_metrics",
      label = "Choose 3 to 8 metrics to peform segmentation :",
      choices = metrics_infos_reactive() %>%
        dplyr::filter(resolution == 20) %>%
        dplyr::filter(type %in% c("lidar", "dendro", "sentinel2")) %>%
        dplyr::pull(name),
      selected = c("z_p95", "z_cv", "z_above2", "sagawi", "fractional_cover_05_2"),
      multiple = TRUE,
      width = "100%",
      options = list(
        plugins = list('remove_button'),
        create = TRUE,
        persist = TRUE,
        maxItems = 8
      )
    )
  })

  # 🟣 Summary metrics list 🟣
  output$summary_metrics_list <- renderUI({
    req(metrics_infos_reactive())
    selectizeInput(
      inputId = "summary_metrics",
      label = "Choose summary metrics :",
      choices = metrics_infos_reactive() %>%
        dplyr::filter(resolution == 20) %>%
        dplyr::filter(type %in% c("lidar", "dendro", "sentinel2")) %>%
        dplyr::pull(name),
      selected = c("dens", "qmdbh", "Vmerch_ha"),
      multiple = TRUE,
      width = "100%",
      options = list(
        'plugins' = list('remove_button'),
        'create' = TRUE,
        'persist' = TRUE
      )
    )
  })

  # 🟣 Mask list 🟣
  output$masks_list <- renderUI({
    req(selected_wd_reactive())

    possibles_masks <- c("roads.shp", "waterbodies.shp", "watercourses.shp")
    shp_lists <- list.files(paste0(selected_wd_reactive(), "/shp"))
    masks_available <- shp_lists[shp_lists %in% possibles_masks]
    masks_available <- sub(".shp", "", masks_available)

    checkboxGroupInput("masks",
                       "Choose the masks to peform segmentation :",
                       choices = setNames(masks_available, str_to_title(masks_available)),
                       selected = c("roads", "waterbodies"))
  })

  # 🟠 Compute enhanced forest resources inventory polygons 🟠
  observeEvent(input$run, {
    # 🟢 Create and set wd 🟢
    segmentation_wd <- paste0(selected_wd_reactive(), "/segmentations/", input$name)
    dir.create(segmentation_wd)

    # 🟢 Best models for imputation 🟢
    read.csv("D:/00_Ontario_eFRI/RMF/results/best_model.csv") -> best_models

    best_models %>%
      pull(knn_vars) %>%
      strsplit(",") %>%
      unlist() %>%
      unique() -> best_models_variables

    # 🟢 Read catalog 🟢
    showNotification("Read catalog", type = "message", duration = 15, session = session)
    st_read(paste0(selected_wd_reactive(), "/ctg/ctg.shp")) -> ctg

    # 🟢 Get epsg from catalog 🟢
    epsg <- st_crs(ctg)

    # 🟢 Extraction circle 🟢
    extraction_circle_reactive() %>%
      st_as_sf() %>%
      st_transform(epsg$wkt) -> extraction_circle

    extraction_circle %>%
      st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
               layer = "extraction_circle",
               quiet = T)

    # 🟢 Save metadata 🟢
    start <- Sys.time()
    cat(c("Start time : ", as.character(start), "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("Segmentation metrics : ", input$segmentation_metrics, "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("Selected summary metrics : ", input$summary_metrics, "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("Coordinates of extraction circle : ", paste0(input$longitude, ", ", input$latitude), "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("Radius of extraction circle (km) : ", input$distance, "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("EPSG : ", epsg$input, "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("Segmentation parameters : ", paste0("thresh = ", input$grm_thresh, " / spec = ", input$grm_spec, " / spat = ", input$grm_spat), "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    cat(c("Segmentation masks : ", paste0(input$masks, collapse = ", "), "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    # 🟢 Read data 🟢
    # Metrics
    showNotification("Read metrics", type = "message", duration = 15, session = session)
    metrics_infos %>%
      filter(name %in% c(best_models_variables,
                         input$segmentation_metrics,
                         input$summary_metrics,
                         "z_p95", "z_above2", "fractional_cover_05_2", "slope", "sagawi")) -> metrics_infos_selected

    metrics_infos_selected %>%
      pull(name) -> metrics_names

    metrics_infos_selected %>%
      pull(path) %>%
      map(rast) %>%
      map(project, epsg$wkt, method = "bilinear") %>%
      map(resample, .[[1]]) %>%
      rast -> metrics

    # Assign metrics names
    names(metrics) <- metrics_names

    # Masks
    if(!is.null(input$masks)){
    showNotification("Read masks", type = "message", duration = 15, session = session)
    paste0(selected_wd_reactive(), "/shp/", input$masks, ".shp") %>%
      map(vect) %>%
      map(project, epsg$wkt) -> masks
    } else {
      showNotification("No masks selected", type = "message", duration = 15, session = session)
      masks <- input$masks
    }

    # Forest inventory polygons (fri)
    showNotification("Read forest inventory polygons", type = "message", duration = 15, session = session)
    st_read(paste0(selected_wd_reactive(), "/shp/PolygonForest.shp")) %>%
      rowid_to_column("id") %>%
      st_transform(epsg$wkt) -> fri_polygons

    # Landcover
    showNotification("Read landcover", type = "message", duration = 15, session = session)
    rast(paste0(selected_wd_reactive(), "/metrics/other/landcover.tif")) %>%
      terra::project(epsg$wkt, method = "near") -> landcover

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
    showNotification("Read forest age", type = "message", duration = 15, session = session)
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_age_2019.tif")) %>%
      terra::project(epsg$wkt, method = "near") -> forest_age_2019

    # Forest_Fire_1985-2020
    showNotification("Read forest fire", type = "message", duration = 15, session = session)
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_fire_1985_2020.tif")) %>%
      terra::project(epsg$wkt, method = "near") -> forest_fire_1985_2020

    # Forest_Harvest_1985-2020
    showNotification("Read forest harvest", type = "message", duration = 15, session = session)
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_harvest_1985_2020.tif")) %>%
      terra::project(epsg$wkt, method = "near") -> forest_harvest_1985_2020


    # 🟢 Clip data 🟢
    showNotification("Clip catalog", type = "message", duration = 15, session = session)
    ctg %<>%
      st_filter(extraction_circle)

    showNotification("Clip metrics", type = "message", duration = 15, session = session)
    metrics %<>%
      crop(extraction_circle) %>%
      mask(extraction_circle)

    if(!is.null(input$masks)){
      showNotification("Clip masks", type = "message", duration = 15, session = session)
      masks %<>%
        map(crop, extraction_circle)
    }

    showNotification("Clip forest inventory polygons", type = "message", duration = 15, session = session)
    fri_polygons %<>%
      st_filter(extraction_circle)

    showNotification("Clip landcover", type = "message", duration = 15, session = session)
    landcover %<>%
      crop(extraction_circle) %>%
      mask(extraction_circle)

    showNotification("Clip forest age", type = "message", duration = 15, session = session)
    forest_age_2019 %<>%
      crop(extraction_circle) %>%
      mask(extraction_circle)

    showNotification("Clip forest fire", type = "message", duration = 15, session = session)
    forest_fire_1985_2020 %<>%
      crop(extraction_circle) %>%
      mask(extraction_circle)

    showNotification("Clip forest harvest", type = "message", duration = 15, session = session)
    forest_harvest_1985_2020 %<>%
      crop(extraction_circle) %>%
      mask(extraction_circle)


    # 🟢 Segmentation 🟢
    showNotification("Perform segmentation", type = "message", duration = 15, session = session)

    eFRI_segmentation(metrics = metrics[[input$segmentation_metrics]],
                      masks = masks,
                      thresh = input$grm_thresh,
                      spec = input$grm_spec,
                      spat = input$grm_spat,
                      method = "bs",
                      clean_nodata = TRUE,
                      output_path = segmentation_wd,
                      output_name = "segmentation",
                      otb_dir = otb_dir) -> segmentation

    segmentation %>%
      st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
               layer = "segmentation",
               quiet = T)

    # 🟢 Imputation 🟢
    showNotification("Peform imputation", type = "message", duration = 15, session = session)

    eFRI_imputation(segmentation = segmentation,
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

    segmentation_imputed %>%
      st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
               layer = "segmentation_imputed",
               quiet = T)

    # 🟢 Build attribute table 🟢
    showNotification("Build attribute table in segmented polygons", type = "message", duration = 15, session = session)

    eFRI_attribute_table(segmentation = segmentation_imputed,
                         metrics = metrics,
                         summary_metrics = paste0(input$summary_metrics, "_median"),
                         landcover = landcover,
                         forest_fire = forest_fire_1985_2020,
                         forest_harvest = forest_harvest_1985_2020,
                         forest_age = forest_age_2019) -> segmentation_data

    segmentation_data %>%
      st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
               layer = "segmentation_data",
               quiet = T)


    # 🟢 Metadata 🟢
    cat(c("Number of forest polygon created : ", nrow(segmentation_data), "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    end <- Sys.time()
    cat(c("End time : ", as.character(end), "----------"),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    elapsed_time <- sub("Time difference of ", "", capture.output(difftime(end, start)))
    cat(c("Elapsed time : ", elapsed_time),
        file = paste0(segmentation_wd, "/metadata.txt"),
        append = TRUE,
        sep = "\n")

    showNotification(paste0("Done after ", elapsed_time, "! A total of ", nrow(segmentation_data), " polygons created for segmentation named ", input$name, "."), type = "message", duration = NULL, session = session)
  })

  # 🟠 Add Segmentation to map 🟠
  observeEvent(input$add_segmentation, {
    segmentation_wd <- paste0(selected_wd_reactive(), "/segmentations/", input$name)
    segmentation <- st_read(dsn = paste0(segmentation_wd, "/data.gpkg"),
                            layer = "segmentation_data",
                            quiet = T) %>%
      st_transform(4326)

    make_popup_generic <- function(x) {
      # x = one row of attributes (data.frame)

      # Format: "colname: value<br>"
      paste0(
        mapply(
          function(name, value) paste0("<b>", name, ":</b> ", value, "<br>"),
          names(x),
          x,
          USE.NAMES = FALSE
        ),
        collapse = ""
      )
    }

    attrs <- sf::st_drop_geometry(segmentation)

    popup_text <- vapply(
      seq_len(nrow(attrs)),
      FUN = function(i) make_popup_generic(attrs[i, ]),
      FUN.VALUE = character(1)
    )

    leafletProxy("map") %>%
      set.view.auto(segmentation) %>%
      addPolygons(data = segmentation,
                  fillOpacity = 0.05,
                  color = "red",
                  weight = 1,
                  popup = popup_text)

  })
}

# Run the application
shinyApp(ui = ui, server = server)
