library(shiny)

# Source function R script
source("~/Capstone/src/dependent_bw.r")
source("~/Capstone/src/plot_shiny.r")


ui <- fillPage(
  shiny::tags$head(
    shiny::tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(30%);
             left: calc(40%);
             font-size: 2em;
             }
             "
      )
    )
  ),
  waiter::use_waiter(),
  titlePanel(h2("Gene Visualization",align="center")),
  # sidebarLayout(
  #   sidebarPanel(
  fillRow(
    flex = c(2,8),
    column(
      width = 12,
      style = "background-color:lightgrey;padding: 20px;margin-left: 10px;",
      textInput("gene", "Enter Gene Name or ID:", value = ""),
      selectInput("run_set", "Select the Run you want to visualize:", 
                  choices = c("mashmap & intronic" = "Run3", "mashmap" = "Run4", "old mapmash & intronic" = "Run2"),
                  selected = 3),
      actionButton("get_gene", "Confirm", class="btn-info btn-block"),
      br(),
      selectInput("dataset", "Select Addional Dataset:", 
                  choices = c("","decoy set" = "d_set", "m1 match set" = "m1_set", "m2 match set" = "m2_set"),
                  selected = NULL),
      br(),
      actionButton("plot","Generate Plot", icon("binoculars"), class = "btn-primary btn-block"),
      br(),
      fluidRow(
        column(6,downloadButton("download_pdf", "PDF", class = "btn-info btn-block")),
        column(6,downloadButton("download_svg", "SVG", class = "btn-block"))
      )
    ),
    # mainPanel(
    column(
      width = 12,
      tabsetPanel(
        tabPanel("Classic Plot",
                 titlePanel(h3(textOutput("title1"),align="center")),
                 plotOutput("classicPlot")),
        tabPanel("Selected Dataset Plot",
                 titlePanel(h3(textOutput("title2"),align="center")),
                 plotOutput("datasetPlot"))
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
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
    req(gene_id)
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
  
  observeEvent(input$plot, {
    req(current()) # ensure tracklist data is available
    
    plot_storage$classic <- function(){
      plotplot(current()[c(1:13,42)],gene_id())
    }
    
    if (input$dataset != "") {
      plot_storage$dataset <- function(){
        subset_data <- switch(input$dataset,
                              "d_set" = current()[c(18:25, 42)],
                              "m1_set" = current()[c(26:33, 42)],
                              "m2_set" = current()[c(34:41, 42)],
                              NULL)
        plotplot(subset_data,gene_id())
      }
    } else {
      plot_storage$dataset <- NULL
    }
  })
  
  output$classicPlot <- renderPlot({
    req(plot_storage$classic)
    plot_storage$classic()
  },res = 80)
  
  output$datasetPlot <- renderPlot({
    req(plot_storage$dataset)
    plot_storage$dataset()
  },res = 100)
  
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
