#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(hivmappr)
library(rstan)
library(rhandsontable)
library(sf)
library(sp)
library(maps)
library(mapproj)
library(dplyr)

library(ggplot2)
library(gridExtra)
library(shiny)
library(leaflet)
library(mapview)


## Load demo data
demo_sh <- sf::as_Spatial(sf::read_sf(system.file("extdata", "mwsh", package="hivmappr")))
demo_csv <- read.csv(system.file("extdata", "mwdf.csv", package="hivmappr"))
demo_mw <- merge(demo_sh, demo_csv)
result_titles <- data.frame(label=c("Prevalence (%)", "ART Coverage (%)", "u_i", "Incidence (per 1000)", "New infections"),
                            var=c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i"),
                            palette=c("Purples", "Blues", "PRGn", "Reds", "Reds"),
                            limit_lower=c(0, 0.3, -0.2, 1e-3, 50),
                            limit_upper=c(0.25, 0.8, 0.2, 8e-3, 1e4),
                            multiplier=c(100, 100, NA, 1000, 1),
                            xintercept=c(0.1, 0.51, NA, 0.0036, NA),
                            logtrans=c(NA, NA, NA, "log", "log10"),
                            rel_min_height=c(0.01, 0.01, 0.01, 0.025, 0.025),
                            stringsAsFactors=FALSE)
transform_func <- function (values, multiplier) {
  if (is.na(multiplier)) return (values)
  round(values*multiplier)
}
th_density <- list(theme_minimal(),
                   theme(plot.title=element_text(hjust=0.5, face="bold"),
                         aspect.ratio=2,
                         axis.text.y=element_blank(),
                         axis.ticks.y = element_line(),
                         axis.ticks.length = unit(2, "pt")))

tempcsvdir <- tempdir()
withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}

shinyServer(function(input, output) {
  values <- reactiveValues(sh=NA, csvdata=NA, mw=NA, data=NA, 
                           imported=file.info(logfile)$mtime,
                           stanlogtext=readLines(logfile)[[1]])
  ## 1. Read in data
  # Read in shapefile
  observeEvent(input$shapefile, {
    if (is.null(input$shapefile)) temp <- 0
    temp_dir <- tempdir()
    # Move uploaded shapefiles to directory
    mv_files <- lapply(1:nrow(input$shapefile), function (i) {
      file.rename(input$shapefile$datapath[i],
                  paste0(temp_dir, "/", input$shapefile$name[i]))
    })
    values$sh <- sf::as_Spatial(sf::read_sf(temp_dir))
    values$mw <- merge(values$sh, values$csvdata)
    showNotification("Your shapefile data have been uploaded.", duration=5)
  })
  observeEvent(input$csv, {
    values$csvdata <- read.csv(input$csv$datapath)
    values$mw <- merge(values$sh, values$csvdata)
    showNotification("Your csv data have been uploaded.", duration=5)
  })
  ## Handsontable
  output$hot <- renderRHandsontable({
    template <- demo_csv[1:2, ]
    template <- data.frame(lapply(template, function (x) {
      output <- switch(class(x), character="", factor="", integer=0, numeric=0)
      rep(output, length(x))
    }))
    rhandsontable(template, readOnly=FALSE)
  })
  observeEvent(input$save, {
    showNotification("Your pasted data table has been saved", duration=5)
    values$csvdata <- hot_to_r(input$hot)
    values$mw <- merge(values$sh, values$csvdata)
  })
  # Demo data
  observeEvent(input$loaddemo, {
    # If the user selects the demo option, then the data in the hivmappr package are used
    showNotification("Demo shapefiles loaded.", duration=5)
    values$mw <- demo_mw
    values$sh <- demo_sh
    values$csvdata <- demo_csv
  })
  output$shapeoutline <- renderPlot({
    if (is.na(values$sh))
      return (print(ggplot()))
    mwpoly <- map_data(values$mw, namefield="district")
    ggplot(mwpoly, aes(long, lat, group = group)) +
      geom_polygon(col="black", fill=NA) +
      coord_map() +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank())
  })
  shapefileplot_leaflet <- reactive({
    if (is.na(values$sh)) return (print(ggplot()))
    var_x <- input$Layers
    geo_vars <- c("district", "region", "zone")
    plot_layer <- paste0('~colorQuantile("Blues", ', var_x,')(', var_x,')')
    discrete_test <- length(unique(values$mw[[var_x]])) < 10
    if (discrete_test) {
      plot_layer <- plot_layer %>% gsub("Quantile", "Factor", .)
    }
    values$mw$mouseover_labels <- paste0(values$mw$district, ":\n", values$mw[[var_x]])
    colorFunc <- colorNumeric
    if (discrete_test) {
      colorFunc <- colorFactor
    }
    leaflet(values$mw) %>% 
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", pal=colorFunc("Blues", domain=values$mw[[var_x]], reverse=TRUE),
                values=eval(parse(text=paste0("~", var_x))),
                title=var_x, opacity=0.8,
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  output$shapefileplot <- renderLeaflet({
    shapefileplot_leaflet()
  })
  shapefileplot_for_download <- reactive({
    shapefileplot_leaflet() %>%
      setView( lng = input$map_center$lng,lat = input$map_center$lat, zoom = input$map_zoom)
  })
  output$dlshapefileplot <- downloadHandler(
    filename = paste0("shapefile_", input$Layers, ".pdf"
    ), 
    content = function(file) {
      mapshot( x = shapefileplot_for_download(), file = file,
               cliprect = "viewport" # the clipping rectangle matches the height & width from the viewing port
               , selfcontained = FALSE # when this was not specified, the function for produced a PDF of two pages: one of the leaflet map, the other a blank page.
      )
    }
  ) 
  ## 4. Fit stan model
  modelfit <- eventReactive(input$modelbutton, {
    values$data <- with(values$mw@data, list(N_reg = length(district),
                                             district = district,
                                             prev_est = prev_survey,
                                             prev_se = prev_survey_se,
                                             anc1_obs = cbind(ancrt_neg = ancrt_n - ancrt_pos,
                                                              ancrt_noart=ancrt_pos - ancrt_art,
                                                              ancrt_art=ancrt_art),
                                             P_i = npos,
                                             R_i = nrecent,
                                             pop15pl_i = pop15pl,
                                             pop15to49_i = pop15to49,
                                             art15pl_i = adultart,
                                             prev_ratio = 1.06,
                                             OmegaT0 = 130 / 365,
                                             sigma_OmegaT = ((142-118)/365)/(2*qnorm(0.975)),
                                             betaT0 = 0.0,
                                             sigma_betaT = 0.00001,
                                             omega = 0.7,
                                             T = 1.0,
                                             Nkappa = 1,
                                             Xkappa = matrix(10, length(district), 1),
                                             sigma_u_sd = 1))
    showNotification("Starting model fitting process... This message will automatically disappear once fitting has completed.",
                     duration=NULL, id="modelfitnotification")
    fit <- sampling(hivmappr:::stanmodels$incidence_rita,
                  data=values$data,
                  control = list(adapt_delta = 0.95))
    est <- summary(fit, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i"))$summary[, "mean"]
    est <- data.frame(param = names(est), value=est)
    est$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", est$param))
    est$district <- factor(values$mw$district[est$district_idx], values$data$district)
    est$region <- factor(values$mw$region[est$district_idx], c("Northern", "Central", "Southern"))
    est$param <- sub("([^\\[]+).*", "\\1", est$param)
    samp <- abind::abind(as.matrix(fit, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i")), along=0)
    names(dimnames(samp)) <- c("model", NA, "param")
    samp <- reshape2::melt(samp)
    samp$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", samp$param))
    samp$district <- factor(values$mw$district[samp$district_idx],
                            rev(values$mw$district))
    samp$region <- factor(values$mw$region[samp$district_idx], c("Northern", "Central", "Southern"))
    samp$param <- sub("([^\\[]+).*", "\\1", samp$param)
    removeNotification("modelfitnotification")
    list(fit=fit, estimates=est, samples=samp)
  })
  output$modelfitprint <- renderText({
    A <- modelfit()
    ""
  })
  output_table <- reactive({
    # Produce table of outputs to view and download
    est <- rstan::summary(modelfit()$fit, c("rho_i", "alpha_i", "lambda_i", "infections_i"))$summary[, c("mean", "2.5%", "97.5%")]
    colnames(est) <- c("mean", "ci_l", "ci_u")
    est <- data.frame(param = rownames(est), est)
    est$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", est$param))
    est$district <- factor(values$mw$district[est$district_idx], values$data$district)
    est$region <- factor(values$mw$region[est$district_idx], c("Northern", "Central", "Southern"))
    est$param <- sub("([^\\[]+).*", "\\1", est$param) %>% sub("\\_i$", "", .)
    est$scale <- c(rho = 100, alpha = 100, lambda = 1000, infections = 1)[est$param]
    est$digits <- c(rho = 1, alpha = 0, lambda = 1, infections = -2)[est$param]
    est$str <- with(est, sprintf(paste0("%.", pmax(digits, 0), "f (%.", pmax(digits, 0), "f, %.", pmax(digits, 0), "f)"),
                                 round(scale*mean, digits),
                                 round(scale*ci_l, digits),
                                 round(scale*ci_u, digits)))
    est$label <- factor(est$param, c("rho", "alpha", "lambda", "infections"),
                        c("Prevalence (%)", "ART coverage (%)", "Incidence (per 1000)", "New infections"))
    reshape2::dcast(est, district+region ~ label, value.var="str")
  })
  output$results_table <- renderTable({
    output_table()
  })
  # Downloadable csv of selected dataset ----
  output$download <- renderUI({
    if (!is.null(output_table())) {
      downloadButton("downloadResults", "Download results")
    }
  })
  output$downloadResults <- downloadHandler(
    filename = function() {
      "results.csv"
    },
    content = function(file) {
      write.csv(output_table(), file, row.names = FALSE)
    }
  )
  # Visualize results
  estplot_leaflet <- reactive({
    estimates <- modelfit()$estimates
    mw_object <- values$mw
    plotting_params <- filter(result_titles, label==input$Results)
    variable <- plotting_params$var
    mw_object@data <- full_join(mw_object@data, estimates[estimates$param==variable, ], by="district")
    plot_layer <- paste0('~colorNumeric("', 
                         plotting_params$palette,'", c(', 
                         plotting_params$limit_lower,', ', 
                         plotting_params$limit_upper,
                         '))(value)')
    mw_object@data$mouseover_labels <- 
      paste0(mw_object@data$district, ":\n", 
             transform_func(mw_object@data$value, plotting_params$multiplier), "%")
    leaflet(mw_object) %>%
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", 
                pal=colorNumeric(plotting_params$palette, 
                                 domain=c(plotting_params$limit_lower, plotting_params$limit_upper), 
                                 reverse=TRUE),
                values=c(plotting_params$limit_lower, plotting_params$limit_upper),
                title=input$Results, opacity=0.8,
                labFormat = labelFormat(transform = function(x) {
                  sort(transform_func(x, plotting_params$multiplier), decreasing = TRUE)
                }))
  })
  output$estplot <- renderLeaflet({
    estplot_leaflet()
  })
  output$densityplot <- renderPlot({
    samples <- modelfit()$samples
    plotting_params <- filter(result_titles, label==input$Results)
    variable <- plotting_params$var
    vertical_line <- geom_vline(color="grey20", linetype="dashed")
    plotdata <- samples[samples$param == variable, ]
    Plot <- ggplot(data=plotdata,
                   aes(x=value, y=district, fill = ..x.. )) +
      ggridges::geom_density_ridges_gradient(rel_min_height=plotting_params$rel_min_height)
    if (!is.na(plotting_params$xintercept)) {
      Plot <- Plot +
        geom_vline(xintercept=plotting_params$xintercept, color="grey20", linetype="dashed")
    }
    Plot <- Plot +
      scale_x_continuous(element_blank(), limits=quantile(plotdata$value, c(0.001, 0.999)),
                         labels=function(x) transform_func(x, plotting_params$multiplier),
                         trans=ifelse(is.na(plotting_params$logtrans), "identity", "log10")) +
      scale_y_discrete(element_blank()) +
      scale_fill_distiller(guide = "none", palette=plotting_params$palette, direction=1) +
      labs(title=input$Results) +
      th_density + theme(axis.text.y=element_text(hjust=1))
    Plot
  }, res=50)
  # Session info
  observe({
    withConsoleRedirect("infoconsole", {
      print(sessionInfo())
    })
  })
})

## TO DO:
# - add download button to Visualize Result panel
# - add download all button
# - work out how to export map
# - display progress bar for R stan or print console output to ui.R
# - add data checks
# - add 'Click here to visualize uploaded data'
# - add MRC logo

# - people upload shapefiles - select variables. on panel 1: add dropdown menu that says pick
# shapefile that corresponds. force the names for datafiles
# - people select other admin levels to visualize
# - when uploading shapefile, populate table automatically
# - populate table when uploading csv file
# - Visualize data - bar plots of prevalence by districts/regions
# 2 versions of models

