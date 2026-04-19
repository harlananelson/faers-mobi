# faers.mobi — Adverse Event Signal Detection (FAERS, precomputed + CSV upload)

box::use(
  shiny[bootstrapPage, moduleServer, NS, tags, div, fluidPage, navbarPage,
        tabPanel, fileInput, reactive, req, observeEvent, h4, p, hr,
        fluidRow, column],
)

box::use(
  app/view/signal_table,
  app/view/signal_timeline,
)

#' @export
ui <- function(id) {
  ns <- NS(id)
  navbarPage(
    title = "faers.mobi — Adverse Event Signal Detection",
    theme = bslib::bs_theme(bootswatch = "flatly"),

    # Primary tab: precomputed signals with caterpillar plot over time
    tabPanel("Signals over time",
      fluidPage(
        div(class = "container-fluid", style = "padding-top: 20px;",
          signal_timeline$ui(ns("timeline"))
        )
      )
    ),

    # Legacy tab: bring-your-own CSV upload
    tabPanel("Upload CSV",
      fluidPage(
        div(class = "container-fluid", style = "padding-top: 20px;",
          fluidRow(
            column(12,
              h4("Bayesian Gamma-Poisson Signal Detection"),
              p("Upload AE report data (CSV with product, event columns) to run ",
                "disproportionality detection live in the browser. ",
                "For a richer view with time-stratified signals across historical ",
                "FAERS data, see the ", tags$em("Signals over time"), " tab."),
              hr()
            )
          ),
          fluidRow(
            column(6,
              fileInput(ns("data_file"), "Upload FAERS/VAERS data (CSV):",
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
          tags$p("Bayesian and frequentist disproportionality analysis over FAERS ",
                 "(FDA Adverse Event Reporting System) data. Signals are precomputed ",
                 "offline via a deterministic R/targets pipeline and served as ",
                 "read-only parquet for interactive exploration."),
          tags$h3("Statistical methods"),
          tags$ul(
            tags$li("GPS/EBGM with 2-component Gamma mixture prior (DuMouchel 1999)"),
            tags$li("PRR + Yates chi-squared (Evans 2001, MHRA criterion)"),
            tags$li("ROR with log-normal CI (van Puijenbroek 2002)"),
            tags$li("BCPNN/IC with 95% credibility bound (Bate 1998, Noren 2006)")
          ),
          tags$h3("Time dimension"),
          tags$p("Per-quarter rolling 4-quarter window with cumulative-fit prior; ",
                 "EWMA smoothing (lambda = 0.3); no Stan hierarchical model in v1."),
          tags$h3("Disclaimer"),
          tags$p("Disproportionate reporting is a statistical pattern, not evidence of ",
                 "causation. Signals are hypotheses requiring further investigation. ",
                 "See ", tags$a(href = "https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html",
                                "FDA FAERS Public Dashboard"), " for official FDA signals.")
        )
      )
    )
  )
}

#' @export
server <- function(id) {
  moduleServer(id, function(input, output, session) {
    signal_timeline$server("timeline")

    uploaded_data <- reactive({
      req(input$data_file)
      utils::read.csv(input$data_file$datapath, stringsAsFactors = FALSE)
    })

    signal_table$server("signals", data = uploaded_data)
  })
}
