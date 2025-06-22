library(shiny)
library(shinydashboard)
library(waiter)
library(shinycssloaders)
library(pryr)
library(profvis)

# Source function R script
source("~/Capstone/src/visual_shiny/utils_data.R")
source("~/Capstone/src/visual_shiny/global_init.R")
source("~/Capstone/src/visual_shiny/plot_tracks.R") 

source("~/Capstone/src/visual_shiny/ui.R")
source("~/Capstone/src/visual_shiny/server.R")

# Run the Shiny app
shinyApp(ui = ui, server = server)

