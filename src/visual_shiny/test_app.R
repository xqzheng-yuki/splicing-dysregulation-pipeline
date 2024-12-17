# Load the Shiny package
library(shiny)
library(log4r)
library(shinydashboard)

my_layout <- function(level, ...) {
  paste(format(Sys.time()), "[", level, "]", ..., sep = " ", collapse = "", "\n")
}
logger <- logger(appenders=console_appender(my_layout))
level(logger) <- "DEBUG"

# Define the User Interface
ui <- dashboardPage(
  dashboardHeader(title="Graph Plot"),
  
  dashboardSidebar(
    numericInput(
      inputId = "num",
      label = "Enter a number:",
      value = 5,
      min = 1,
      max = 100
    ),
    sidebarSearchForm("gene","get_gene",label = "Enter Name or ID",icon = shiny::icon("ok",lib = "glyphicon")),
    actionButton("plot", "Plot"),
    textInput("run_set",NULL,value="Run3"),
    helpText("The graph will plot values from 1 to the input number and their squares.")),
  
  dashboardBody(
    actionButton("prev", "Previous"),
    actionButton("nex", "Next"),
    plotOutput("plot1"),
    plotOutput("plot2"),
    plotOutput("classicPlot")
  )
)

# Define the Server Logic
server <- function(input, output) {
  state <- reactiveValues(results = list(), current_index = 0)
  plot_storage <- reactiveValues(classic = NULL, dataset = NULL)
  
  reactive_data <- reactive({
    x <- 1:input$num
    y <- x^2
    list(x = x, y = y)
  })
  
  gene_id <- eventReactive(input$get_gene, {
    req(input$gene)  # Ensure input$gene is available
    # Check if the gene input is already an ID or needs to be converted from a name
    if (is_gene_name(input$gene)) {
      get_gene_id(input$gene)  # If input is a gene name, use get_gene_id() to convert
    } else {
      input$gene  # If input is already a gene ID, return it as-is
    }
  })
  current <- eventReactive(input$plot, {
    req(gene_id())
    req(input$run_set)
    id <- showNotification("Getting track data...", duration = NULL, closeButton = FALSE,type="warning")
    on.exit(removeNotification(id), add = TRUE)
    waiter::Waiter$new(id = "classicPlot")$show()
    tracklist(gene_id(),input$run_set)
  })
  add_new_plots <- function(state, classic_plot_func, dataset_plot_func) {
    state$results <- append(state$results, list(
      list(classic = classic_plot_func, dataset = dataset_plot_func)
    ))
    state$current_index <- length(state$results) # Set to the latest added
    info(logger,paste0("The state is ",length(state$results)))
    state
  }
  
  observeEvent(input$plot, {
    req(reactive_data())
    req(current())
    data <- reactive_data()
    # plot_storage$classic <- function(){
    #   plot(data$x, data$y, type = "l", col = "red", lwd = 2,
    #        xlab = "x", ylab = "y",
    #        main = "Line Plot of x vs y = x^2")}
    plot_storage$classic <- function(){
      plotplot(current()[c(1:13,42)],gene_id())
    }
    plot_storage$dataset <- function(){
      plot(data$x, data$y, type = "p", pch = 19, col = "blue",
           xlab = "x", ylab = "y",
           main = "Points Plot of x vs y = x^2")}
    state <- add_new_plots(state, plot_storage$classic, plot_storage$dataset)
    state$current_index <- length(state$results)
  })
  
  observeEvent(input$prev, {
    if (state$current_index > 1) {
      state$current_index <- state$current_index - 1
      # Update plot_storage with the selected set
      plot_storage$classic <- state$results[[state$current_index]]$classic
      plot_storage$dataset <- state$results[[state$current_index]]$dataset
      info(logger,paste0("Now it is on graph ",state$current_index))
    }
  })
  observeEvent(input$nex, {
    if (state$current_index < length(state$results)) {
      state$current_index <- state$current_index + 1
      # Update plot_storage with the selected set
      plot_storage$classic <- state$results[[state$current_index]]$classic
      plot_storage$dataset <- state$results[[state$current_index]]$dataset
      info(logger,paste0("Now it is on graph ",state$current_index))
    }
  })
  
  output$plot1 <- renderPlot({
    req(plot_storage$classic)
    plot_storage$classic()
  },res = 80)
  
  output$plot2 <- renderPlot({
    req(plot_storage$dataset)
    plot_storage$dataset()
  },res = 80)
  
}

# Run the Shiny App
shinyApp(ui = ui, server = server)
