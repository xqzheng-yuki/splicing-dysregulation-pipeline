library(shiny)
library(shinydashboard)
library(waiter)
library(shinycssloaders)
library(pryr)
library(profvis)

# Source function R script
source("~/Capstone/src/visual_shiny/dependent_bw.r")
source("~/Capstone/src/visual_shiny/plot_shiny.r")

titlecontent <- readLines(file.path("www", "guide.txt"), warn = FALSE)

ui <- dashboardPage(
  dashboardHeader(title=span(
    "Gene Visualization",
    span(
      `data-toggle` = "tooltip", `data-placement` = "right",`data-html` = "true",
      title = HTML(titlecontent),
      icon("info-circle")
    )
  )),
  dashboardSidebar(
      # textInput("gene", "Enter Gene Name or ID:", value = ""),
      sidebarSearchForm("gene","get_gene",label = "Enter Name or ID",icon = shiny::icon("ok",lib = "glyphicon")),
      helpText("Please clilck the check mark after input your interest gene."),
      selectInput("run_set", "Select the Run you want to visualize:", 
                  choices = c("mashmap & intronic" = "Run3", "mashmap" = "Run4", "old mapmash & intronic" = "Run2"),
                  selected = "Run3"),
      selectInput("dataset", "Select Addional Dataset:", 
                  choices = c("decoy set" = "d_set", "m1 match set" = "m1_set", "m2 match set" = "m2_set"),
                  selected = "d_set"),
      br(),
      actionButton("plot","Generate Plot", icon("binoculars"), width = '200px', class = "btn-primary btn-block"),
      br(),
      fillRow(
        column(12,downloadButton("download_pdf", "PDF", class = "btn-info btn-block")),
        column(12,downloadButton("download_svg", "SVG", class = "btn-block"))
      )
    ),
  dashboardBody(
    shiny::tags$head(
      shiny::tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      shiny::tags$script(src = "scripts.js")
    ),
    waiter::use_waiter(),
    shiny::tags$div(
      class = "right-align-btn",
      actionButton("prev", "Previous"),
      actionButton("nex", "Next")
    ),
      tabBox(width = 12,
        tabPanel("Classic Plot",
                 titlePanel(h3(textOutput("title1"),align="center")),
                 plotOutput("classicPlot",height = 900)),
        tabPanel("Selected Dataset Plot",
                 titlePanel(h3(textOutput("title2"),align="center")),
                 plotOutput("datasetPlot"))
      )
    )
  )

# Define Server
server <- function(input, output, session) {
  
  # Reactive state to store historical plots
  state <- reactiveValues(results = list(), current_index = 0)
  # Reactive values to store plots
  plot_storage <- reactiveValues(classic = NULL, dataset = NULL)
  
  # Determine if input is a gene name or gene ID
  gene_id <- eventReactive(input$get_gene, {
    req(input$gene)  # Ensure input$gene is available
    # Check if the gene input is already an ID or needs to be converted from a name
    if (is_gene_name(input$gene)) {
      get_gene_id(input$gene)  # If input is a gene name, use get_gene_id() to convert
    } else {
      input$gene  # If input is already a gene ID, return it as-is
    }
  })
  
  # For get title for plot
  output$title1 <- output$title2 <- renderText({
    req(gene_id())
    gene_id()
  })

  # Generate the tracklist based on the gene ID
  current <- eventReactive(input$plot, {
    req(gene_id())
    req(input$run_set)
    id <- showNotification("Getting track data...", duration = NULL, closeButton = FALSE,type="warning")
    on.exit(removeNotification(id), add = TRUE)
    waiter::Waiter$new(id = "classicPlot")$show()
    tracklist(gene_id(),input$run_set)
  })
  
  # Generate plots
  observeEvent(input$plot, {
    req(current()) # ensure tracklist data is available
    print(paste("Object size of current track:", object.size(current)))
    # Plot for main
    plot_storage$classic <- function(){
      plotplot(current()[c(1:13,42)],gene_id())
    }
    
    # Plot for additional
    plot_storage$dataset <- function(){
      subset_data <- switch(input$dataset,
                            "d_set" = current()[c(18:25, 42)],
                            "m1_set" = current()[c(26:33, 42)],
                            "m2_set" = current()[c(34:41, 42)])
      plotplot(subset_data,gene_id())
    }
    
    # Add the new set of plots to the state for history
    state <- add_new_plots(state, plot_storage$classic, plot_storage$dataset)
    # Update the current index to the latest set
    state$current_index <- length(state$results)
    print(paste("Number of plots stored:", length(state$results)))
    print(paste("Memory usage after adding plots:", pryr::mem_used()))
    print(paste("Object size of result:", object.size(state$results)))
  })
  
  # Navigating-Up
  observeEvent(input$prev, {
    showNotification("You click previous plot...",duration = 0.5, closeButton = TRUE,type="warning")
    if (state$current_index > 1) {
      state$current_index <- state$current_index - 1
      # Update plot_storage with the selected set
      plot_storage$classic <- state$results[[state$current_index]]$classic %||% function() NULL
      plot_storage$dataset <- state$results[[state$current_index]]$dataset %||% function() NULL
      info(logger,paste0("Now it is on graph ",state$current_index))
    }
  })
  
  # Navigating-Down
  observeEvent(input$nex, {
    showNotification("You click for next plot...",duration = 0.5, closeButton = TRUE,type="warning")
    if (state$current_index < length(state$results)) {
      state$current_index <- state$current_index + 1
      # Update plot_storage with the selected set
      plot_storage$classic <- state$results[[state$current_index]]$classic %||% function() NULL
      plot_storage$dataset <- state$results[[state$current_index]]$dataset %||% function() NULL
      info(logger,paste0("Now it is on graph ",state$current_index))
    }
  })
  
  # Render plot-1
  output$classicPlot <- renderPlot({
    req(plot_storage$classic)
    print(paste("Memory before rendering classicPlot:", pryr::mem_used()))
    plot_storage$classic()
    print(paste("Memory after rendering classicPlot:", pryr::mem_used()))
  },res = 80)
  
  # Render plot-2
  output$datasetPlot <- renderPlot({
    req(plot_storage$dataset)
    print(paste("Memory before rendering datasetPlot:", pryr::mem_used()))
    plot_storage$dataset()
    print(paste("Memory after rendering datasetPlot:", pryr::mem_used()))
  },res = 100)
  
  # Download-pdf
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste("plots-", gene_id(), Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      if (!is.null(plot_storage$classic)) {
        plot_storage$classic()
      }
      if (!is.null(plot_storage$dataset)) {
        plot_storage$dataset()
      }
      dev.off()
    }
  )
  
  #Download-plot
  output$download_svg <- downloadHandler(
    filename = function() {
      paste("plots-", gene_id(), Sys.Date(), ".svg", sep = "")
    },
    content = function(file) {
      svg(file)
      if (!is.null(plot_storage$classic)) {
        plot_storage$classic()
      }
      if (!is.null(plot_storage$dataset)) {
        plot_storage$dataset()
      }
      dev.off()
    }
  )
}

# Run the Shiny app

  shinyApp(ui = ui, server = server)

