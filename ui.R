library(shiny)
library(readxl)
library(tidyverse)
library(plotly)
library(limma)
library(plyr)

shinyUI(fluidPage(
  titlePanel(title=h1("MLPA interactive plotter")),
  sidebarLayout(
    sidebarPanel(position="left",
                 helpText("This application will fit linear models on gene expression data from multiplex 
                          ligation-dependent probe amplification experiments (see Eldering 2003)
                          with repeated measurements and visualize log fold changes between two
                          conditions in a plotly volcano plot.", align = "center"),
                 helpText(        "You can either", align = "center"),
                 fileInput("fo", "Upload your dataset in CSV format"),
                 helpText("OR", align = "center"),
                 tags$br(),
                 actionButton("settempl", "Load the template dataset",align="center"),
                 tags$br(),      
                 "The template file with instructions can be found",a("here",href="template.xlsx"),
                 tags$hr(),
                 uiOutput("dropdown1"),
                 uiOutput("dropdown2"),
                 htmlOutput("instruction")
                 ),
    mainPanel(
      plotlyOutput("volcplot"))
  )
)
  
)