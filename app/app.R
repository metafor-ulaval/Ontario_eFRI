# 🟡🟡 Load libraries🟡🟡 ----
library(shiny)
library(shinyjs)
library(shinyBS)
library(leaflet)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(magrittr)
library(sf)
library(terra)
library(mapedit)
library(callr)
# pak::pak("metafor-ulaval/eFRItools")
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

# 🌐 Memory given to OTB and GDAL, based on the memory available on the computer 🌐
memory_budget <- function(fraction = 0.6) {
  avail_mb <- ps::ps_system_memory()$avail / 1024^2
  ram_mb <- max(floor(avail_mb * fraction), 1024)

  list(ram_mb = ram_mb,
       gdal_cache_mb = floor(ram_mb / 4))
}

# 🌐 Find OTB, first in the app folder, then in the working directory 🌐
find_otb_dir <- function(app_dir, wd) {
  c(paste0(app_dir, "/softwares/OTB-9.1.0-Win64/bin"),
    paste0(wd, "/softwares/OTB-9.1.0-Win64/bin")) %>%
    keep(dir.exists) %>%
    head(1)
}





# 🟡🟡 Parameters 🟡🟡 ----
# The app folder is the folder of this file (shiny::runApp sets it as working directory)
app_dir <- normalizePath(getwd(), winslash = "/")

# When the app is inside the deliverable folder, the deliverable folder is the default working directory
default_wd <- normalizePath(paste0(app_dir, "/.."), winslash = "/")
if (!all(dir.exists(paste0(default_wd, c("/metrics", "/shapefiles"))))) default_wd <- NULL

# When launched with the launcher, closing the browser stops the app
desktop_mode <- Sys.getenv("EFRI_DESKTOP") == "1"

n_steps <- 7

# One log file per launch of the app, in app/logs
log_dir <- paste0(app_dir, "/logs")
dir.create(log_dir, showWarnings = FALSE)
app_log_file <- paste0(log_dir, "/eFRI_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log")

write_log <- function(...) {
  cat(paste0(format(Sys.time(), "%H:%M:%S"), " ", paste0(...), "\n"), file = app_log_file, append = TRUE)
}

# Run an action, writing any error and its call stack in the log instead of stopping the app
run_logged <- function(action_name, session, expr, on_error = NULL) {
  tryCatch(
    withCallingHandlers(expr, error = function(e) {
      if (inherits(e, "shiny.silent.error")) return() # req() stops silently, it is not an error
      write_log("ERROR in ", action_name, " : ", conditionMessage(e), "\nCall stack :\n",
                paste0("  ", vapply(sys.calls(), function(x) paste(deparse(x, nlines = 1), collapse = ""), character(1)), collapse = "\n"))
    }),
    error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      if (!is.null(on_error)) on_error()
      showNotification(paste0("Unexpected error in ", action_name, " : ", conditionMessage(e), " See the log in ", app_log_file, "."),
                       type = "error", duration = NULL, session = session)
      NULL
    })
}

write_log("App started | R ", R.version$major, ".", R.version$minor, " | app folder : ", app_dir, " | working directory : ", if (is.null(default_wd)) "none" else default_wd)





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
                            # Working directory
                            actionButton("wd_browse", "Choose a working directory"),
                            verbatimTextOutput("selected_wd"),
                            uiOutput("forest_list"),
                            # Parameters
                            uiOutput("segmentation_metrics_list"), # Peut-être mettre un min et un max a cocher
                            uiOutput("summary_metrics_list"), # Peut-être mettre un min et un max a cocher
                            uiOutput("masks_list"),
                            strong("Choose generic region merging parameters :"),
                            numericInput("grm_thresh", "Threshold (0-100)", 50, min = 0, max = 100, step = 1),
                            # strong("Maximum heterogeneity allowed when merging, controlling how fine or coarse the segmentation is"),
                            # "Small threshold → fine segmentation (many small regions).",
                            # "Large threshold → coarse segmentation (fewer, larger regions).",
                            numericInput("grm_spec", "Weight of spectral homogeneity (0-1)", 0.5, min = 0, max = 1, step = 0.1),
                            # strong("Weight given to spectral similarity (pixel values, band means) when deciding whether regions should merge."),
                            # "If high → segmentation will mostly respect spectral values.",
                            # "If low → spectral similarity matters little, other criteria (shape) dominate.",
                            numericInput("grm_spat", "Weight of spatial homogeneity (0-1)", 0.5, min = 0, max = 1, step = 0.1),
                            # strong("Weight given to shape similarity (compactness and smoothness) when deciding whether regions should merge."),
                            # "If high → the algorithm favors compact, smooth regions even if spectral similarity is weaker.",
                            # "If low → region boundaries will mostly follow spectral homogeneity.",
                            # Action button to show plots and compute statistics
                            textInput("name", "Choose output name :", ""),
                            actionButton("run", "Compute enhanced forest resources inventory polygons"),
                            hidden(actionButton("cancel", "Cancel computation", class = "btn-danger")),
                            hidden(verbatimTextOutput("run_log")),
                            actionButton("add_segmentation", "Add enhanced forest resources inventory polygons to map")))
      )
    ),

    # 🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵

    # 🔵 Interface principale 🔵
    shiny::mainPanel(
        div(class = "sidebar-panel",
            # Parameters
            editModUI(id = "map_draw_module", height = "90vh", width = "80vw"))
    )

    # 🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵🔵

  )
)





# 🟡🟡 Server 🟡🟡 ----
server <- function(input, output, session) {

  # 🟢 = Steps
  # 🟠 = Reactive components
  # 🟣 = Action button components

  # 🟠 Working directory via the native OS dialog 🟠
  wd_value <- reactiveVal(default_wd)

  observeEvent(input$wd_browse, {
    start_dir <- if (is.null(wd_value())) "" else wd_value()

    path <- tryCatch({
      if (.Platform$OS.type == "windows") {
        utils::choose.dir(default = start_dir,
                          caption = "Select the eFRI working directory")
      } else {
        tcltk::tk_choose.dir(default = start_dir,
                             caption = "Select the eFRI working directory")
      }
    }, error = function(e) NA_character_)

    # NA = user cancelled
    if (length(path) == 1 && !is.na(path)) {
      wd_value(normalizePath(path, winslash = "/", mustWork = FALSE))
    }
  })

  selected_wd_reactive <- reactive({
    req(wd_value())
    wd_value()
  })

  # Directory displayed in text
  output$selected_wd <- renderText({
    if (is.null(wd_value())) "No directory selected" else wd_value()
  })

  # 🟢 Map initialisation 🟢
  select_area <- callModule(mapedit::editMod,
                            id = "map_draw_module",
                            leafmap = leaflet() %>% addProviderTiles("Esri.WorldImagery", group = "Satellite"),
                            sf = TRUE,
                            record = FALSE)

  # 🟣 Forest list 🟣
  output$forest_list <- renderUI({
    req(selected_wd_reactive())
    radioButtons(
      inputId = "forest",
      label = "Choose the forest where to peform segmentation :",
      choices = list.files(paste0(selected_wd_reactive(), "/metrics")),
      selected = character(0),
      width = "100%",
    )
  })

  # 🟠 Reactive ctg 🟠
  ctg_reactive <- reactive({
    req(input$forest, selected_wd_reactive())

    st_read(paste0(selected_wd_reactive(), "/shapefiles/", input$forest, "/ctg.shp"), quiet = TRUE)
  })

  # 🟠 Update ctg 🟠
  observe({
    req(ctg_reactive())

    ctg_reactive()%>%
      st_buffer(1) %>%
      st_union() %>%
      st_transform(4326) -> ctg_clean

    leafletProxy("map_draw_module-map") %>%
      set_view_auto(ctg_clean) %>%
      clearGroup("ctg") %>%
      addPolygons(data = ctg_clean,
                  group = "ctg",
                  color = "white",
                  weight = 2,
                  fillOpacity = 0)

  })

  # 🟠 Read metrics informations 🟠
  metrics_infos_reactive <- reactive({
    req(input$forest, selected_wd_reactive())
    list.files(paste0(selected_wd_reactive(), "/metrics/", input$forest), full.names = TRUE) %>%
      read_metrics()
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
      selected = c("z_p95", "z_above2", "z_skew"),
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
      selected = c("vmerch_ha", "qmdbh", "dens", "ba_ha"),
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
    req(input$forest, selected_wd_reactive())

    possibles_masks <- c("roads.shp", "waterbodies.shp")
    shp_lists <- list.files(paste0(selected_wd_reactive(), "/shapefiles/", input$forest))
    masks_available <- shp_lists[shp_lists %in% possibles_masks]
    masks_available <- sub(".shp", "", masks_available)

    checkboxGroupInput("masks",
                       "Choose the masks to peform segmentation :",
                       choices = setNames(masks_available, str_to_title(masks_available)),
                       selected = c("roads", "waterbodies"))
  })

  # 🟠 Background computation state 🟠
  run_state <- reactiveValues(process = NULL,
                              progress = NULL,
                              segmentation_wd = NULL,
                              name = NULL)

  stop_run <- function() {
    if (!is.null(run_state$progress)) run_state$progress$close()
    run_state$process <- NULL
    run_state$progress <- NULL
    shinyjs::enable("run")
    shinyjs::hide("cancel")
    shinyjs::hide("run_log")
  }

  # 🟠 Compute enhanced forest resources inventory polygons 🟠
  observeEvent(input$run, run_logged("Compute", session, on_error = function() {
    # Remove the output folder if nothing was written in it, and reset the buttons
    segmentation_wd <- paste0(isolate(selected_wd_reactive()), "/segmentation/", isolate(input$forest), "/automated_", isolate(input$name))
    if (dir.exists(segmentation_wd) && length(list.files(segmentation_wd)) == 0) unlink(segmentation_wd, recursive = TRUE)
    if (!is.null(isolate(run_state$process)) && isolate(run_state$process)$is_alive()) isolate(run_state$process)$kill_tree()
    stop_run()
  }, {

    write_log("Compute clicked | forest : ", input$forest, " | name : ", input$name,
              " | segmentation metrics : ", paste(input$segmentation_metrics, collapse = ", "),
              " | summary metrics : ", paste(input$summary_metrics, collapse = ", "),
              " | masks : ", paste(input$masks, collapse = ", "),
              " | thresh / spec / spat : ", input$grm_thresh, " / ", input$grm_spec, " / ", input$grm_spat,
              " | subset area drawn : ", !is.null(select_area()$finished))

    # 🟢 Condition 0 🟢
    if(!is.null(run_state$process)){
      showNotification("A computation is already running.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 1 🟢
    if(is.null(wd_value())){
      showNotification("Choose working directory.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 2 🟢
    if(is.null(input$forest)){
      showNotification("Choose forest.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 3 🟢
    if(!isTruthy(input$name)){
      showNotification("Choose output name.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 3b 🟢
    if(!grepl("^[A-Za-z0-9_-]+$", input$name)){
      showNotification("Output name can only contain letters (without accents), numbers, '-' and '_'.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Set wd 🟢
    paste0(selected_wd_reactive(), "/segmentation/", input$forest, "/automated_", input$name) -> segmentation_wd

    # 🟢 Condition 4 🟢
    if(file.exists(paste0(segmentation_wd, "/data.gpkg"))){
      showNotification("This segmentation already exist.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 5 🟢
    if(length(input$segmentation_metrics) < 3){
      showNotification("Choose at least 3 segmentation metrics.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 5b 🟢
    if(!isTruthy(input$grm_thresh) || input$grm_thresh <= 0 || input$grm_thresh > 100){
      showNotification("Threshold must be greater than 0 and lower or equal to 100.", type = "message", duration = 15, session = session)
      return()
    }

    if(!isTruthy(input$grm_spec) || input$grm_spec <= 0 || input$grm_spec > 1){
      showNotification("Weight of spectral homogeneity must be greater than 0 and lower or equal to 1.", type = "message", duration = 15, session = session)
      return()
    }

    if(!isTruthy(input$grm_spat) || input$grm_spat <= 0 || input$grm_spat > 1){
      showNotification("Weight of spatial homogeneity must be greater than 0 and lower or equal to 1.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 6 🟢
    otb_dir <- find_otb_dir(app_dir, selected_wd_reactive())
    if(length(otb_dir) == 0){
      showNotification("Orfeo ToolBox (OTB-9.1.0-Win64) was not found in the softwares folder.", type = "error", duration = 15, session = session)
      return()
    }

    # 🟢 Get epsg from dem 🟢
    metrics_infos_reactive() %>%
      filter(name == "dem") %>%
      pull(path) %>%
      rast() %>%
      st_crs() -> epsg

    # 🟢 Condition 7 🟢
    extraction_area <- NULL
    if(!is.null(select_area()$finished)){
      select_area()$finished %>%
        st_as_sf() %>%
        st_transform(epsg) -> extraction_area

      ctg_reactive() %>%
        st_transform(epsg) -> ctg

      if(!any(st_intersects(ctg, extraction_area, sparse = FALSE))){
        showNotification("Subset area is outside of catalog, please place area within catalog or remove it.", type = "message", duration = 15, session = session)
        return()
      }
    } else {
      showNotification("No subset area selected, the whole forest management unit will be processed. This can take a long time.", type = "warning", duration = 30, session = session)
    }

    if(!is.null(extraction_area)){
      write_log("Subset area | features : ", nrow(extraction_area), " | bbox : ", paste(round(sf::st_bbox(extraction_area)), collapse = ", "),
                " | columns : ", paste(names(extraction_area), collapse = ", "))
    }

    # 🟢 Create wd 🟢
    dir.create(segmentation_wd, recursive = TRUE)
    write_log("Output folder created : ", segmentation_wd)

    # 🟢 Launch computation in background 🟢
    memory <- memory_budget()
    write_log("Memory (MB) | OTB : ", memory$ram_mb, " | GDAL cache : ", memory$gdal_cache_mb)

    params <- list(wd = selected_wd_reactive(),
                   forest = input$forest,
                   segmentation_wd = segmentation_wd,
                   segmentation_metrics = input$segmentation_metrics,
                   summary_metrics = input$summary_metrics,
                   masks = input$masks,
                   thresh = input$grm_thresh,
                   spec = input$grm_spec,
                   spat = input$grm_spat,
                   extraction_area = extraction_area,
                   otb_dir = otb_dir,
                   ram_mb = memory$ram_mb,
                   gdal_cache_mb = memory$gdal_cache_mb)

    run_state$process <- callr::r_bg(function(app_dir, p) {
                                       source(paste0(app_dir, "/R/pipeline.R"))
                                       run_efri_pipeline(p)
                                     },
                                     args = list(app_dir = app_dir, p = params),
                                     stdout = paste0(segmentation_wd, "/log.txt"),
                                     stderr = "2>&1",
                                     supervise = TRUE)
    write_log("Background computation started | pid : ", run_state$process$get_pid())

    run_state$segmentation_wd <- segmentation_wd
    run_state$name <- input$name
    run_state$progress <- shiny::Progress$new(session, min = 0, max = n_steps)
    run_state$progress$set(value = 0, message = "Starting computation")

    shinyjs::disable("run")
    shinyjs::show("cancel")
    shinyjs::show("run_log")
  }))

  # 🟠 Follow background computation 🟠
  observe(run_logged("Follow computation", session, on_error = stop_run, {
    req(run_state$process)
    invalidateLater(1000)

    # Last step written in the log file of the segmentation (==== Step 1/7 : Read metrics ====)
    log_file <- paste0(run_state$segmentation_wd, "/log.txt")
    if(file.exists(log_file)){
      steps <- grep("^==== Step [0-9]+/[0-9]+ : .* ====$", readLines(log_file, warn = FALSE), value = TRUE)
      if(length(steps) != 0){
        step <- as.numeric(sub("^==== Step ([0-9]+)/.*$", "\\1", tail(steps, 1)))
        run_state$progress$set(value = step - 0.5,
                               message = gsub("^==== | ====$", "", tail(steps, 1)))
      }
    }

    if(run_state$process$is_alive()) return()

    result <- tryCatch(run_state$process$get_result(), error = function(e) e)

    if(inherits(result, "error")){
      message <- if (!is.null(result$parent)) conditionMessage(result$parent) else conditionMessage(result)
      cat(c("Failed : ", message),
          file = paste0(run_state$segmentation_wd, "/metadata.txt"),
          append = TRUE,
          sep = "\n")
      write_log("Computation failed : ", message)
      showNotification(paste0("The computation failed : ", message, " See log.txt in the output folder for details."),
                       type = "error", duration = NULL, session = session)
    } else {
      write_log("Computation done | polygons : ", result$n_polygons, " | elapsed time : ", result$elapsed_time)
      showNotification(paste0("Done after ", result$elapsed_time, ". A total of ", result$n_polygons, " polygons created for segmentation named ", run_state$name, "."),
                       type = "message", duration = NULL, session = session)
    }

    stop_run()
  }))

  # 🟠 Last lines of the log 🟠
  output$run_log <- renderText({
    req(run_state$process)
    invalidateLater(1000)

    log_file <- paste0(run_state$segmentation_wd, "/log.txt")
    if(!file.exists(log_file)) return("")
    paste(tail(readLines(log_file, warn = FALSE), 8), collapse = "\n")
  })

  # 🟣 Cancel computation 🟣
  observeEvent(input$cancel, {
    req(run_state$process)

    run_state$process$kill_tree()
    unlink(run_state$segmentation_wd, recursive = TRUE)
    write_log("Computation cancelled : ", run_state$segmentation_wd)
    showNotification(paste0("Computation of segmentation named ", run_state$name, " cancelled."), type = "message", duration = 15, session = session)

    stop_run()
  })

  # 🟠 Stop computation (and app in desktop mode) when the browser is closed 🟠
  session$onSessionEnded(function() {
    process <- isolate(run_state$process)
    if (!is.null(process) && process$is_alive()) process$kill_tree()
    if (desktop_mode) stopApp()
  })

  # 🟠 Add Segmentation to map 🟠
  observeEvent(input$add_segmentation, {

    # 🟢 Condition 1 🟢
    if(is.null(wd_value())){
      showNotification("Choose working directory.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 2 🟢
    if(is.null(input$forest)){
      showNotification("Choose forest.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 3 🟢
    if(!isTruthy(input$name)){
      showNotification("Choose output name.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Set wd 🟢
    paste0(selected_wd_reactive(), "/segmentation/", input$forest, "/automated_", input$name) -> segmentation_wd

    # 🟢 Condition 4 🟢
    if(!dir.exists(segmentation_wd)){
      showNotification("This segmentation doesn't exist.", type = "message", duration = 15, session = session)
      return()
    }

    # 🟢 Condition 5 🟢
    metadata_file <- paste0(segmentation_wd, "/metadata.txt")
    metadata <- if (file.exists(metadata_file)) readLines(metadata_file, warn = FALSE) else character(0)
    failed <- which(metadata == "Failed : ")

    if(length(failed) != 0){
      showNotification(paste0("This segmentation failed : ", paste(metadata[-seq_len(failed[1])], collapse = " ")),
                       type = "error", duration = NULL, session = session)
      return()
    }

    # 🟢 Condition 6 🟢
    layer_names <- if (file.exists(paste0(segmentation_wd, "/data.gpkg"))) st_layers(paste0(segmentation_wd, "/data.gpkg"))$name else character(0)

    if(!any(layer_names == "segmentation_data_imputed")){
      showNotification("This segmentation is not completed.", type = "message", duration = 15, session = session)
      return()
    }

    st_read(dsn = paste0(segmentation_wd, "/data.gpkg"),
                            layer = "segmentation_data_imputed",
                            quiet = T) %>%
      st_transform(4326) -> segmentation

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

    leafletProxy("map_draw_module-map") %>%
      set_view_auto(segmentation) %>%
      clearGroup("segmentation") %>%
      addPolygons(data = segmentation,
                  group = "segmentation",
                  fillOpacity = 0.05,
                  color = "red",
                  weight = 1,
                  popup = popup_text)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
