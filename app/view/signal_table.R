# Signal results table module

box::use(
  shiny[moduleServer, NS, reactive, req, renderTable, tableOutput,
        selectInput, sliderInput, actionButton, observeEvent,
        tagList, tags, div, fluidRow, column, h3, p, hr],
)

#' @export
ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4,
        selectInput(ns("product_filter"), "Filter by product:",
                    choices = NULL, multiple = TRUE)
      ),
      column(4,
        sliderInput(ns("threshold"), "EB05 signal threshold:",
                    min = 0.5, max = 5, value = 1, step = 0.25)
      ),
      column(4,
        div(style = "margin-top: 25px;",
          actionButton(ns("run"), "Detect Signals", class = "btn-primary btn-lg")
        )
      )
    ),
    hr(),
    h3("Signal Detection Results"),
    tableOutput(ns("results_table"))
  )
}

#' @export
server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    box::use(
      app/logic/signal_engine[run_signal_detection, get_signals],
    )

    results <- reactive({
      req(input$run)
      d <- data()

      if (!is.null(input$product_filter) && length(input$product_filter) > 0) {
        d <- d[d$product %in% input$product_filter, ]
      }

      run_signal_detection(d, threshold = input$threshold)
    }) |> shiny::bindEvent(input$run)

    output$results_table <- renderTable({
      req(results())
      signals <- get_signals(results())
      if (nrow(signals) == 0) {
        return(data.frame(Message = "No signals detected at this threshold."))
      }
      signals$rr <- round(signals$rr, 2)
      signals$eb05 <- round(signals$eb05, 3)
      signals$eb50 <- round(signals$eb50, 3)
      signals$eb95 <- round(signals$eb95, 3)
      signals
    })
  })
}
