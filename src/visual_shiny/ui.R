titlecontent <- readLines(file.path("~/Capstone/src/visual_shiny/www", "guide.txt"), warn = FALSE)

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
    sidebarSearchForm("gene","get_gene",label = "Enter Name or ID",icon = shiny::icon("ok",lib = "glyphicon")),
    p(strong("Hit'ENTER' or the check mark"),style="text-align:center"),
    p(strong("after entering gene."),style="text-align:center"),
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
                    plotOutput("datasetPlot",height = 600))
    )
  )
)