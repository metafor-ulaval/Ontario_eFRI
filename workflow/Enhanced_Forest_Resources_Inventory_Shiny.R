# 🟡🟡 Load libraries🟡🟡 ----
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyFiles)
library(leaflet)
library(leaflet.extras)
library(tidyverse)
library(magrittr)
library(sf)
library(terra)
#library(knitr)
library(exactextractr)
library(eFRItools)





# 🟡🟡 Paramters 🟡🟡 ----
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
                            uiOutput("metric_list"), # Peut-être mettre un min et un max a cocher
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

  # 🟢 Read metrics informations 🟢
  list.files(paste0(selected_wd_reactive(), "/metrics"), full.names = TRUE) %>%
    read_metrics() -> metrics_infos

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
    # Catalog
    st_read(paste0(selected_wd_reactive(), "/ctg/ctg.shp")) %>%
      st_transform(input$epsg) -> ctg

    # Subset metrics in a 2 km radius circle in the centroid of area ################### CLIP
    ctg %>% ################### CLIP
      st_bbox() %>% ################### CLIP
      st_as_sfc() %>% ################### CLIP
      st_centroid() %>% ################### CLIP
      st_buffer(1000) %>% ################### CLIP
      st_as_sf() -> subset_circle ################### CLIP

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
      map(project, input$epsg, method = "bilinear") %>%
      map(resample, .[[1]]) %>%
      rast -> metrics

    # Assign metrics names
    names(metrics) <- metrics_names

    ################### CLIP
    metrics %<>% ################### CLIP
      crop(subset_circle) %>% ################### CLIP
      mask(subset_circle) ################### CLIP

    # Masks
    paste0(selected_wd_reactive(), "/shp/", input$masks, ".shp") %>%
      map(vect) %>%
      map(project, input$epsg) -> masks

    # Forest inventory polygons (fri)
    st_read(paste0(selected_wd_reactive(), "/shp/PolygonForest.shp")) %>%
      rowid_to_column("id") %>%
      st_transform(input$epsg) -> fri_polygons

    ################### CLIP
    fri_polygons %<>% ################### CLIP
      st_filter(subset_circle) ################### CLIP

    # Best models for imputation
    read.csv(paste0(selected_wd_reactive(), "/results/best_model.csv")) %>%
      pull(knn_vars) %>%
      strsplit(",") %>%
      unlist() %>%
      unique() -> best_models_variables

    # Landcover
    rast(paste0(selected_wd_reactive(), "/metrics/other/landcover.tif")) -> landcover

    # Forest_Age_1985-2020
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_age_2019.tif")) -> forest_age_2019

    # Forest_Fire_1985-2020
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_fire_1985_2020.tif")) -> forest_fire_1985_2020

    # Forest_Harvest_1985-2020
    rast(paste0(selected_wd_reactive(), "/metrics/other/forest_harvest_1985_2020.tif")) -> forest_harvest_1985_2020




    # 🟢 Segmentation 🟢
    showNotification("Perform segmentation", type = "message", duration = 15, session = session)

    eFRI_segmentation(metrics = metrics,
                      masks = masks,
                      thresh = input$grm_thresh,
                      spec = input$grm_spec,
                      spat = input$grm_spat,
                      method = "bs",
                      output_path = selected_wd_reactive(),
                      output_name = input$name,
                      otb_dir = otb_dir)

    # 🟢 Imputation 🟢
    showNotification("Peform imputation", type = "message", duration = 15, session = session)

    eFRI_imputation(segmentation = st_read(paste0(selected_wd_reactive(), input$name)),
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

    eFRI_attribute_table(segmentation = segmentation,
                         metrics = metrics,
                         landcover = landcover,
                         forest_fire = forest_fire_1985_2020,
                         forest_harvest = forest_harvest_1985_2020,
                         forest_age = forest_age_2019) -> segmentation_data

    showNotification("Done !", type = "message", duration = 15, session = session)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
