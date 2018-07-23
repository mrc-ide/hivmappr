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

th_density <- list(theme_minimal(),
                   theme(plot.title=element_text(hjust=0.5, face="bold"),
                         aspect.ratio=2,
                         axis.text.y=element_blank(),
                         axis.ticks.y = element_line(),
                         axis.ticks.length = unit(2, "pt")))

tempcsvdir <- tempdir()

shinyServer(function(input, output) {
  values <- reactiveValues(sh=NA, csvdata=NA, mw=NA, data=NA)
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
    values$mw$mouseover_labels <- paste0(values$mw$district, ":\n", values$mw[[var_x]])
    leaflet(values$mw) %>% 
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", pal=colorNumeric("Blues", domain=values$mw[[var_x]], reverse=TRUE),
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
  estprevplot_leaflet <- reactive({
    estimates <- modelfit()$estimates
    mw_object <- values$mw
    mw_object@data <- full_join(mw_object@data, estimates[estimates$param == "rho_i",], by="district")
    plot_layer <- '~colorNumeric("Purples", c(0, 0.25))(value)'
    mw_object@data$mouseover_labels <- 
      paste0(mw_object@data$district, ":\n", round(mw_object@data$value*100), "%")
    leaflet(mw_object) %>%
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", pal=colorNumeric("Purples", domain=c(0, 0.25), reverse=TRUE),
                values=c(0, 0.25),
                title="Prevalence (%)", opacity=0.8,
                labFormat = labelFormat(transform = function(x) sort(round(100*x), decreasing = TRUE)))
  })
  output$estprevplot <- renderLeaflet({
    estprevplot_leaflet()
  })
  estartcoverageplot_leaflet <- reactive({
    estimates <- modelfit()$estimates
    mw_object <- values$mw
    mw_object@data <- full_join(mw_object@data, estimates[estimates$param=="alpha_i",], by="district")
    plot_layer <- '~colorNumeric("Blues", c(0.3, 0.8))(value)'
    mw_object@data$mouseover_labels <- 
      paste0(mw_object@data$district, ":\n", round(mw_object@data$value*100), "%")
    leaflet(mw_object) %>%
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", pal=colorNumeric("Blues", domain=c(0.3, 0.8), reverse=TRUE),
                values=c(0.3, 0.8),
                title="ART Coverage (%)", opacity=0.8,
                labFormat = labelFormat(transform = function(x) sort(round(100*x), decreasing = TRUE)))
  })
  output$estartcoverageplot <- renderLeaflet({
    estartcoverageplot_leaflet()
  })
  estuiplot_leaflet <- reactive({
    estimates <- modelfit()$estimates
    mw_object <- values$mw
    mw_object@data <- full_join(mw_object@data, estimates[estimates$param=="u_i",], by="district")
    plot_layer <- '~colorNumeric("PRGn", c(-0.2, 0.2))(value)'
    mw_object@data$mouseover_labels <- 
      paste0(mw_object@data$district, ":\n", mw_object@data$value)
    leaflet(mw_object) %>%
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", pal=colorNumeric("PRGn", domain=c(-0.2, 0.2), reverse=TRUE),
                values=c(-0.2, 0.2),
                title="u_i", opacity=0.8,
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  output$estuiplot <- renderLeaflet({
    estuiplot_leaflet()
  })
  estincidenceplot_leaflet <- reactive({
    estimates <- modelfit()$estimates
    mw_object <- values$mw
    mw_object@data <- full_join(mw_object@data, estimates[estimates$param=="lambda_i",], by="district")
    plot_layer <- '~colorNumeric("Reds", c(1, 8)/1e3)(value)'
    mw_object@data$mouseover_labels <- 
      paste0(mw_object@data$district, ":\n", round(1000*mw_object@data$value))
    leaflet(mw_object) %>%
      addPolygons(color="white", opacity=0.5, weight=1,
                  fillColor=eval(parse(text=plot_layer)), fillOpacity=1.0, 
                  highlightOptions=highlightOptions(color="blue", bringToFront = TRUE, sendToBack=TRUE),
                  label=~mouseover_labels) %>%
      addLegend("topright", pal=colorNumeric("Reds", domain=c(1, 8)/1e3, reverse=TRUE),
                values=c(1, 8)/1e3,
                title="Incidence / 1000", opacity=0.8,
                labFormat = labelFormat(transform = function(x) sort(round(1000*x), decreasing = TRUE)))
  })
  output$estincidenceplot <- renderLeaflet({
    estincidenceplot_leaflet()
  })
  output$densityplots <- renderPlot({
    samples <- modelfit()$samples
    panA <- ggplot(data=samples[samples$param == "rho_i", ],
                   aes(x=value, y=district, fill = ..x.. )) +
      ggridges::geom_density_ridges_gradient(rel_min_height=0.01) +
      geom_vline(xintercept=0.10, color="grey20", linetype="dashed") +
      scale_x_continuous(element_blank(), labels=function(x) round(100*x, 1)) +
      scale_y_discrete(element_blank()) +
      scale_fill_distiller(guide = "none", palette="Purples", direction=1, trans="log10") +
      labs(title="Prevalence (%)") +
      th_density + theme(axis.text.y=element_text(hjust=1))
    panB <- ggplot(data=samples[samples$param == "alpha_i", ],
                   aes(x=value, y=district, fill = ..x.. )) +
      ggridges::geom_density_ridges_gradient(rel_min_height=0.01) +
      geom_vline(xintercept=0.51, color="grey20", linetype="dashed") +
      scale_x_continuous(element_blank(), labels=function(x) round(100*x, 1)) +
      scale_y_discrete(element_blank()) +
      scale_fill_distiller(guide = "none", palette="Blues", direction=1) +
      labs(title="ART Coverage (%)") +
      th_density
    panC <- ggplot(data=samples[samples$param == "lambda_i", ],
                   aes(x=value, y=district, fill = ..x.. )) +
      ggridges::geom_density_ridges_gradient(rel_min_height=0.025) +
      geom_vline(xintercept=0.0036, color="grey20", linetype="dashed") +
      scale_x_log10(element_blank(),
                    breaks=c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02),
                    limits=c(0.0006, 0.025), labels=function(x) round(1000*x, 1)) +
      scale_y_discrete(element_blank()) +
      scale_fill_distiller(guide = "none", palette="Reds", direction=1) +
      labs(title="Incidence / 1000 (log)") +
      th_density
    panD <- ggplot(data=samples[samples$param == "infections_i", ],
                   aes(x=value, y=district, fill = ..x.. )) +
      ggridges::geom_density_ridges_gradient(rel_min_height=0.025) +
      scale_x_log10(element_blank(), limits=c(50, 10000), labels=function(x) round(x)) +
      scale_y_discrete(element_blank()) +
      scale_fill_distiller(guide = "none", palette="Reds", direction=1) +
      labs(title="New infections", fontface="bold") +
      th_density
    grid.arrange(panA, panB, panC, panD, ncol=4, widths=c(1.35, 1, 1, 1))
  }, res=50)
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
  output$results <- renderTable({
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
  output$sessioninfo <- renderPrint(capture.output(sessionInfo()))
})


