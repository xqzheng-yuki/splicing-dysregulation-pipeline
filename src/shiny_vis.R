library(shiny)
library(Gviz)

# Source the provided R scripts (assuming they define functions such as plot_gene and plot_supplymental)
source("~/Capstone/src/visualization_bw.r")
source("~/Capstone/src/dependent_bw.r")

# Define UI for the Shiny app
ui <- fluidPage(
  titlePanel("Gene Visualization"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene", "Gene to Visualize:", value = "Enter gene name here", placeholder = "Enter gene name here"),
      actionButton("plotButton", "Generate Plot"),
      downloadButton("downloadPlot", "Download Plot as PDF")
    ),
    mainPanel(
      plotOutput("genePlot"),
      plotOutput("supplementPlot")
    )
  )
)

# Define server logic for the Shiny app
server <- function(input, output, session) {
  # Reactive expression to generate the plots when the button is clicked
  plot_data <- eventReactive(input$plotButton, {
    req(input$gene)
    
    # Create temporary files for saving the plots
    gene_plot_file <- tempfile(fileext = ".pdf")
    supplement_plot_file <- tempfile(fileext = ".pdf")
    
    # Call your custom plotting functions from the sourced scripts
    # Save the plots to temporary files instead of directly to disk
    pdf(gene_plot_file)
    plot_gene(input$gene)
    dev.off()
    
    pdf(supplement_plot_file)
    plot_supplymental(input$gene)
    dev.off()
    
    list(gene_plot_file = gene_plot_file, supplement_plot_file = supplement_plot_file)
  })
  
  # Render the gene plot in the UI
  output$genePlot <- renderPlot({
    req(plot_data())
    pdf_file <- plot_data()$gene_plot_file
    plot_content <- pdf_text(pdf_file)
    plot(plot_content)
  })
  
  # Render the supplement plot in the UI
  output$supplementPlot <- renderPlot({
    req(plot_data())
    pdf_file <- plot_data()$supplement_plot_file
    plot_content <- pdf_text(pdf_file)
    plot(plot_content)
  })
  
  # Allow users to download the plots as a PDF
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("gene_plot_", input$gene, ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      gene_pdf <- plot_data()$gene_plot_file
      supplement_pdf <- plot_data()$supplement_plot_file
      plot(gene_pdf)
      plot(supplement_pdf)
      dev.off()
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
