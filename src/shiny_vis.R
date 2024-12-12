library(shiny)
library(shinydashboard)
library(waiter)
library(shinycssloaders)

# Source function R script
source("~/Capstone/src/dependent_bw.r")
source("~/Capstone/src/plot_shiny.r")
css <- '
        .shiny-notification {
             position:fixed;
             top: calc(30%);
             left: calc(40%);
             font-size: 2em;
        }
        .main-content{
            overflow-y: auto;
            height: 100vh
        }
        .tooltip {
            pointer-events: none;
        }
        .tooltip > .tooltip-inner {
              pointer-events: none;
              background-color: #E0E0E0;
              color: black;
              padding: 10px;
              font-size: 15px;
              font-style: normal;
              text-align: left;
              margin-left: 0;
              max-width: 300px;
        }
        .tooltip > .arrow::before {
              border-right-color: black;
        }
'
js <- "
$(function () {
  $('[data-toggle=tooltip]').tooltip()
})
"
titlecontent <- "Usage Instruction:<br>
    <strong>1. Enter Name or ID:</strong> Input the desired name or ID into the text box at the top and confirm by clicking the checkmark icon.<br>
    <strong>2. Select Run:</strong> Use the dropdown to choose the run you want to visualize.<br>
    <strong>3. Select Additional Dataset (Optional):</strong> Use the second dropdown to select an additional dataset, if needed.<br>
    <strong>4. Generate Plot:</strong> Click the blue 'Generate Plot' button to create the visualization.<br>
    <strong>5. Download Output:</strong> Click the 'PDF' or 'SVG' buttons to download the generated plot in your preferred format.
  "
ui <- dashboardPage(
  dashboardHeader(title=span(
    "Gene Visualization",
    span(
      `data-toggle` = "tooltip", `data-placement` = "right",`data-html` = "true",
      title = titlecontent,
      icon("info-circle")
    )
  )),
  # titlePanel(h2("Gene Visualization",align="center")),
  # sidebarLayout(
  #   sidebarPanel(
  # fillRow(
  #   flex = c(2,8),
  #   column(
  #     width = 12,
  #     style = "background-color:lightgrey;padding: 20px;margin-left: 10px;",
  dashboardSidebar(
      # textInput("gene", "Enter Gene Name or ID:", value = ""),
      sidebarSearchForm("gene","get_gene",label = "Enter Name or ID",icon = shiny::icon("ok",lib = "glyphicon")),
      selectInput("run_set", "Select the Run you want to visualize:", 
                  choices = c("mashmap & intronic" = "Run3", "mashmap" = "Run4", "old mapmash & intronic" = "Run2"),
                  selected = "Run3"),
      # column(12,actionButton("get_gene", "Confirm", class="btn-info btn-block")),
      selectInput("dataset", "Select Addional Dataset:", 
                  choices = c("","decoy set" = "d_set", "m1 match set" = "m1_set", "m2 match set" = "m2_set"),
                  selected = NULL),
      br(),
      actionButton("plot","Generate Plot", icon("binoculars"), width = '200px', class = "btn-primary btn-block"),
      br(),
      fillRow(
        column(12,downloadButton("download_pdf", "PDF", class = "btn-info btn-block")),
        column(12,downloadButton("download_svg", "SVG", class = "btn-block"))
      )
    ),
    # mainPanel(
    # column(
    #   width = 12,
  dashboardBody(
    shiny::tags$head(
      shiny::tags$style(HTML(css)),
      tags$script(HTML(js))
    ),
    waiter::use_waiter(),
      tabBox(width = 12,
        tabPanel("Classic Plot",
                 titlePanel(h3(textOutput("title1"),align="center")),
                 plotOutput("classicPlot",height = 900)),
        tabPanel("Selected Dataset Plot",
                 titlePanel(h3(textOutput("title2"),align="center")),
                 shinycssloaders::withSpinner(plotOutput("datasetPlot")))
      )
    )
  )
# )

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
