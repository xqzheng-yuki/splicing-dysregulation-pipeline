library(shiny)

ui <- fluidPage(
  titlePanel("Dynamic UI for File Selection"),
  sidebarLayout(
    sidebarPanel(
      actionButton("refresh", "Refresh File List"),  # To refresh file list
      uiOutput("file_selector")  # Dynamic UI for file selection
    ),
    mainPanel(
      textOutput("selected_file")
    )
  )
)

server <- function(input, output, session) {
  # Reactive expression to list files in the directory
  file_list <- reactive({
    input$refresh  # This makes it re-run when the button is clicked
    list.dirs(path = "/mnt/gtklab01/xiaoqing", full.names = FALSE,recursive = FALSE)  # Update with your directory
  })
  
  # Dynamic UI output
  output$file_selector <- renderUI({
    files <- file_list()
    if (length(files) == 0) {
      return("No files found in the directory.")
    }
    selectInput("selected_file", "Choose a file:", choices = files)
  })
  
  # Display selected file
  output$selected_file <- renderText({
    req(input$selected_file)
    paste("You selected:", input$selected_file)
  })
}

shinyApp(ui, server)
