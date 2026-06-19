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
set_view_auto <- function(map,
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
                            # uiOutput("masks_list"),
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
                            actionButton("run", "Compute enhanced forest resources inventory polygons")))
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

  # 🟠 Reactive analysis radius 🟠
  extraction_circle_reactive <- reactive({
    req(input$longitude, input$latitude, input$distance)
    sf::st_sfc(sf::st_point(c(input$longitude, input$latitude)), crs = 4326) %>%
      sf::st_buffer(input$distance*1000) -> extraction_circle

    extraction_circle
  })

  # 🟢 Map initialisation 🟢
  output$map <- renderLeaflet({
    leaflet::leaflet() %>%
      leaflet::addTiles() %>%
      leaflet::addProviderTiles(providers$Esri.WorldImagery)
  })

  # 🟠 Update metrics area and zoom the circle at it 🟠
  observe({
    req(metrics_area_reactive())

    metrics_area_reactive() %>%
      sf::st_transform(4326) -> metrics_area_4326

    leaflet::leafletProxy("map") %>%
      set_view_auto(metrics_area_4326) %>%
      leaflet::addPolygons(data = metrics_area_4326,
                           color = "white",
                           weight = 2,
                           fillOpacity = 0)

    metrics_area_4326 %>%
      sf::st_centroid() %>%
      sf::st_coordinates() -> metrics_area_coordinates
    updateNumericInput(session, "longitude", value = metrics_area_coordinates[1])
    updateNumericInput(session, "latitude", value = metrics_area_coordinates[2])
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

    leaflet::leafletProxy("map") %>%
      leaflet::clearGroup("recherche") %>%
      set_view_auto(extraction_circle_reactive()) %>%
      leaflet::addCircleMarkers(lng = input$longitude,
                                lat = input$latitude,
                                group = "recherche",
                                color = "turquoise",
                                radius = 0.5,
                                weight = 0.5,
                                fillOpacity = 1) %>%
      leaflet::addPolygons(data = extraction_circle_reactive(),
                           color = "turquoise",
                           weight = 2,
                           fillOpacity = 0,
                           group = "recherche",
                           options = pathOptions(interactive = FALSE))
  })

  # 🟠 Read metrics informations 🟠
  metrics_infos_reactive <- reactive({
    req(selected_wd_reactive())

    list.files(selected_wd_reactive(), full.names = T) %>%
      eFRItools::read_metrics() -> metrics_infos

    metrics_infos
  })

  # 🟠 Reactive area 🟠
  metrics_area_reactive <- reactive({
    req(metrics_infos_reactive())

    metrics_infos_reactive() %>%
      dplyr::slice(1) %>%
      dplyr::pull(path) %>%
      terra::rast() -> metric

    metric[!is.na(metric)] <- 1

    metric %>%
      terra::as.polygons() %>%
      sf::st_as_sf() %>%
      nngeo::st_remove_holes() -> metrics_area

    metrics_area
  })

  # 🟣 Metric list 🟣
  output$segmentation_metrics_list <- renderUI({
    req(metrics_infos_reactive())
    selectizeInput(
      inputId = "segmentation_metrics",
      label = "Choose 3 to 8 metrics to peform segmentation :",
      choices = metrics_infos_reactive() %>%
        dplyr::pull(name),
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

  # # 🟣 Mask list 🟣
  # output$masks_list <- renderUI({
  #   req(selected_wd_reactive())
  #
  #   possibles_masks <- c("roads.shp", "waterbodies.shp", "watercourses.shp")
  #   shp_lists <- list.files(paste0(selected_wd_reactive(), "/shp"))
  #   masks_available <- shp_lists[shp_lists %in% possibles_masks]
  #   masks_available <- sub(".shp", "", masks_available)
  #
  #   checkboxGroupInput("masks",
  #                      "Choose the masks to peform segmentation :",
  #                      choices = setNames(masks_available, str_to_title(masks_available)),
  #                      selected = c("roads", "waterbodies"))
  # })

  # 🟠 Compute enhanced forest resources inventory polygons 🟠
  observeEvent(input$run, {
    req(metrics_infos_reactive(), metrics_area_reactive())

    # 🟢 Create and set wd 🟢
    segmentation_wd <- paste0(selected_wd_reactive(), "/", input$name)
    dir.create(segmentation_wd)

    # 🟢 Get epsg from metrics area 🟢
    epsg <- sf::st_crs(metrics_area_reactive())

    # 🟢 Extraction circle 🟢
    extraction_circle_reactive() %>%
      sf::st_as_sf() %>%
      sf::st_transform(epsg$wkt) -> extraction_circle

    extraction_circle %>%
      sf::st_write(dsn = paste0(segmentation_wd, "/data.gpkg"),
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

    # cat(c("Segmentation masks : ", paste0(input$masks, collapse = ", "), "----------"),
    #     file = paste0(segmentation_wd, "/metadata.txt"),
    #     append = TRUE,
    #     sep = "\n")

    # 🟢 Read data 🟢
    # Metrics
    showNotification("Read metrics", type = "message", duration = 15, session = session)
    metrics_infos %>%
      dplyr::filter(name %in% c(input$segmentation_metrics)) -> metrics_infos_selected

    metrics_infos_selected %>%
      dplyr::pull(name) -> metrics_names

    metrics_infos_selected %>%
      dplyr::pull(path) %>%
      purrr::map(terra::rast) %>%
      purrr::map(terra::project, epsg$wkt, method = "bilinear") %>%
      purrr::map(terra::resample, .[[1]]) %>%
      terra::rast() -> metrics

    # Assign metrics names
    names(metrics) <- metrics_names

    # # Masks
    # if(!is.null(input$masks)){
    # showNotification("Read masks", type = "message", duration = 15, session = session)
    # paste0(selected_wd_reactive(), "/shp/", input$masks, ".shp") %>%
    #   map(vect) %>%
    #   map(project, epsg$wkt) -> masks
    # } else {
    #   showNotification("No masks selected", type = "message", duration = 15, session = session)
    #   masks <- input$masks
    # }

    # 🟢 Clip data 🟢
    showNotification("Clip metrics", type = "message", duration = 15, session = session)
    metrics %<>%
      terra::crop(extraction_circle) %>%
      terra::mask(extraction_circle)

    # if(!is.null(input$masks)){
    #   showNotification("Clip masks", type = "message", duration = 15, session = session)
    #   masks %<>%
    #     map(crop, extraction_circle)
    # }

    # 🟢 Segmentation 🟢
    showNotification("Perform segmentation", type = "message", duration = 15, session = session)

    eFRI_segmentation(metrics = metrics[[input$segmentation_metrics]],
                      # masks = masks,
                      masks = NULL,
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

    # 🟢 Metadata 🟢
    cat(c("Number of forest polygon created : ", nrow(segmentation), "----------"),
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

    showNotification(paste0("Done after ", elapsed_time, "! A total of ", nrow(segmentation), " polygons created for segmentation named ", input$name, "."), type = "message", duration = NULL, session = session)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
