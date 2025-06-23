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
    plot_storage$classic <- {
      plotplot(current()[c(1,2,5:12,41)],gene_id())
      recordPlot()
    }
    # plot_storage$classic
    # Plot for additional
    plot_storage$dataset <- {
      subset_data <- switch(input$dataset,
                            "d_set" = current()[c(1,2,17:24, 41)],
                            "m1_set" = current()[c(1,2,25:32, 41)],
                            "m2_set" = current()[c(1,2,33:40, 41)])
      plotplot(subset_data,gene_id())
      recordPlot()
    }
    # plot_storage$dataset
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
      plot_storage$classic <- state$results[[state$current_index]]$classic
      plot_storage$dataset <- state$results[[state$current_index]]$dataset
      info(logger,paste0("Now it is on graph ",state$current_index))
    }
  })
  
  # Navigating-Down
  observeEvent(input$nex, {
    showNotification("You click for next plot...",duration = 0.5, closeButton = TRUE,type="warning")
    if (state$current_index < length(state$results)) {
      state$current_index <- state$current_index + 1
      # Update plot_storage with the selected set
      plot_storage$classic <- state$results[[state$current_index]]$classic
      plot_storage$dataset <- state$results[[state$current_index]]$dataset
      info(logger,paste0("Now it is on graph ",state$current_index))
    }
  })
  
  # Render plot-1
  output$classicPlot <- renderPlot({
    req(plot_storage$classic)
    replayPlot(plot_storage$classic)
  },res = 90)
  
  # Render plot-2
  output$datasetPlot <- renderPlot({
    req(plot_storage$dataset)
    replayPlot(plot_storage$dataset)
  },res = 90)
  
  # Download-pdf
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste("plots-", gene_id(), Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      if (!is.null(plot_storage$classic)) {
        replayPlot(plot_storage$classic)
      }
      if (!is.null(plot_storage$dataset)) {
        replayPlot(plot_storage$dataset)
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
        replayPlot(plot_storage$classic)
      }
      if (!is.null(plot_storage$dataset)) {
        replayPlot(plot_storage$dataset)
      }
      dev.off()
    }
  )
}