# faers.mobi — Vaccine Adverse Event Signal Detection (VAERS)

box::use(
  shiny[bootstrapPage, moduleServer, NS, tags, div, fluidPage, navbarPage,
        tabPanel, fileInput, reactive, req, observeEvent, h4, p, hr,
        fluidRow, column],
)

box::use(
  app/view/signal_table,
)

#' @export
ui <- function(id) {
  ns <- NS(id)
  navbarPage(
    title = "faers.mobi — Vaccine Safety Signal Detection",
    theme = bslib::bs_theme(bootswatch = "flatly"),
    tabPanel("Signal Detection",
      fluidPage(
        div(class = "container-fluid", style = "padding-top: 20px;",
          fluidRow(
            column(12,
              h4("Bayesian Gamma-Poisson Signal Detection for Vaccines"),
              p("Upload VAERS data (CSV with vaccine, adverse event columns) to detect ",
                "disproportionate reporting signals using a 2-component Gamma mixture posterior."),
              hr()
            )
          ),
          fluidRow(
            column(6,
              fileInput(ns("data_file"), "Upload VAERS data (CSV):",
                        accept = ".csv")
            ),
            column(6,
              p(tags$strong("Expected columns:"), " product, event"),
              p("Each row = one adverse event report.")
            )
          ),
          signal_table$ui(ns("signals"))
        )
      )
    ),
    tabPanel("About",
      fluidPage(
        div(style = "max-width: 800px; margin: 40px auto;",
          tags$h2("About faers.mobi"),
          tags$p("Bayesian disproportionality analysis for vaccine adverse event reporting ",
                 "using data from the Vaccine Adverse Event Reporting System (VAERS)."),
          tags$h3("Statistical Method"),
          tags$ul(
            tags$li("Prior: 2-component Gamma mixture fitted via EM algorithm"),
            tags$li("Posterior: full Gamma mixture (not EBGM approximation)"),
            tags$li("Signal: EB05 (5th percentile of posterior) > threshold"),
            tags$li("Reference: DuMouchel (1999), American Statistician 53(3)")
          )
        )
      )
    )
  )
}

#' @export
server <- function(id) {
  moduleServer(id, function(input, output, session) {
    uploaded_data <- reactive({
      req(input$data_file)
      utils::read.csv(input$data_file$datapath, stringsAsFactors = FALSE)
    })

    signal_table$server("signals", data = uploaded_data)
  })
}
