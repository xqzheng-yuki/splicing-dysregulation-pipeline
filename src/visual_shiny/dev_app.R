library(shiny)
library(shinydashboard)
library(waiter)
library(shinycssloaders)
library(pryr)
library(profvis)

# Source function R script
source("~/Capstone/src/visual_shiny/dependent_bw.r")
source("~/Capstone/src/visual_shiny/plot_shiny.r")
source("~/Capstone/src/visual_shiny/major_plot_function.r") 

source("~/Capstone/src/visual_shiny/dev_ui.R")
source("~/Capstone/src/visual_shiny/dev_server.R")

# Run the Shiny app
shinyApp(ui = ui, server = server)

