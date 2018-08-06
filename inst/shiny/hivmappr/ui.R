#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(rhandsontable)
library(leaflet)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  includeCSS("hivmappr.css"),
  titlePanel("hivmappr"),
  navlistPanel(
    tabPanel("0. Home",
             p("If you are interested in testing this app, press the button below to load demo data. You can skip steps 1 and 2 listed in the panel on the left."),
             p("Otherwise, upload your own shapefiles (panel 1) and csv file containing the data (panel 2)."),
             actionButton("loaddemo", "Load demo files"),
             fluidRow(column(6, leafletOutput("initshapefileplotdemo")),
                      column(6, tableOutput("initcsvfile")))
             ),
    # Upload file
    tabPanel("1. Upload shape files",
             fileInput("shapefile", "Upload shapefiles", multiple = TRUE),
             fluidRow(column(6, leafletOutput("initshapefileplot")))
             ),
    tabPanel("2. Input subnational data",
             p("You can either (A) copy and paste data into the spreadsheet below, or (B) you have the option of uploading a csv file containing the data."),
             strong("Option (A)"),
             tags$ol(
               tags$li("Copy table from excel. Do not include the column headings."), 
               tags$li("Go to web browser containing the hivmappr shiny app."), 
               tags$li("Select the first cell in the spreadsheet."),
               tags$li("Paste data by pressing ctrl+v or cmd+v.")
             ),
             actionButton("save", "Save data for next steps"),
             rHandsontableOutput("hot"),
             br(),
             strong("Option (B)"),
             fileInput("csv", "Upload csv file", multiple = FALSE,
                       accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
             ),
    tabPanel("3. Visualize data", 
             fluidRow(column(10, leafletOutput("shapefileplot"), 
                             downloadButton("dlshapefileplot", "Download map"),
                             plotOutput("visualizedatabars"),
                             downloadButton("dlbars", "Download plot")),
                      column(2, selectInput("Layers", "Choose a layer:",
                                            list("pop15pl", "pop15to49", "prev_survey", "prev_survey_se", 
                                                 "nsamp", "npos", "nrecent", "adultart",
                                                 "anc_clients", "ancrt_n", "ancrt_pos", "ancrt_art"))))
             ),
    tabPanel("4. Fit model",
             p("Click the button below to fit the model to your dataset. This may take a few minutes."),
             p("After the model fitting process has completed, you will see your data below and have the option to download it as a csv file."),
             p("Visualizations of the results will appear in Panel 5."),
             actionButton("modelbutton", "Run model"),
             br(),
             textOutput("modelfitprint"), 
             uiOutput("download"),
             tableOutput("results_table")
            ),
    tabPanel("5. Visualize results",
             selectInput("Results", "Show estimates of:", 
                                         as.list(c("Prevalence (%)", "ART Coverage (%)", 
                                                   "u_i", "Incidence (per 1000)", 
                                                   "New infections"))),
             fluidRow(column(7, leafletOutput("estplot")),
                      column(5, plotOutput("densityplot", width="90%"))
                      ),
             fluidRow(column(7, downloadButton("downloadestplot", "Download map")),
                      column(5, downloadButton("downloaddensityplot", "Download estimate densities"))
             )),
    tabPanel("sessionID()", pre(id="infoconsole")),
    widths=c(2,8)
  )
))