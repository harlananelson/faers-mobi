# Signal Timeline — caterpillar plot of credible intervals over time for a
# single (drug, event) pair.
#
# Reads precomputed signals from /srv/shiny-server/faers-mobi/data/signals.parquet
# (produced offline by signal-compute). Never calls safetysignal live — all
# signals are precomputed on the local GPU box.

box::use(
  shiny[NS, moduleServer, tagList, selectizeInput, plotOutput, renderPlot,
        req, reactive, tags, div, h4, p, hr, fluidRow, column, withProgress,
        updateSelectizeInput],
  arrow[open_dataset],
  dplyr[filter, collect, pull, arrange, distinct, `%>%`, mutate, case_when],
)

SIGNALS_PATH <- "data/signals.parquet"
DRUG_DICT_PATH <- "data/drug_dictionary.parquet"
EVENT_DICT_PATH <- "data/event_dictionary.parquet"

#' @export
ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        h4("Signal credible intervals over time"),
        p(
          "Pick a drug and an adverse event. The caterpillar plot shows the ",
          "95% credible interval (EB05–EB95) of the Bayesian Gamma-Poisson ",
          "relative reporting rate per quarter, with the EWMA-smoothed point ",
          "estimate overlaid. Threshold lines mark EBGM = 1 (null) and ",
          "EBGM = 2 (FDA signal criterion)."
        ),
        hr()
      )
    ),
    fluidRow(
      column(6, selectizeInput(ns("drug"), "Drug (RxNorm):",
                                choices = NULL, multiple = FALSE)),
      column(6, selectizeInput(ns("event"), "Event (MedDRA PT):",
                                choices = NULL, multiple = FALSE))
    ),
    fluidRow(column(12, plotOutput(ns("timeline"), height = "500px"))),
    fluidRow(
      column(12,
        tags$small(tags$em(
          "Disclaimer: disproportionate reporting is a statistical pattern, ",
          "not evidence of causation. Signals are hypotheses requiring further ",
          "investigation."
        ))
      )
    )
  )
}

#' @export
server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Lazy-load signals dataset. `arrow::open_dataset` is memory-mapped; does
    # not materialize the whole 200-500 MB parquet into R memory.
    signals <- reactive({
      if (!file.exists(SIGNALS_PATH)) return(NULL)
      open_dataset(SIGNALS_PATH)
    })

    # Populate drug/event selects from dictionaries
    shiny::observe({
      if (file.exists(DRUG_DICT_PATH)) {
        drug_choices <- open_dataset(DRUG_DICT_PATH) %>%
          distinct(.data$rxnorm_name) %>%
          arrange(.data$rxnorm_name) %>%
          collect() %>%
          pull(.data$rxnorm_name)
        updateSelectizeInput(session, "drug", choices = drug_choices,
                             selected = "rofecoxib", server = TRUE)
      }
      if (file.exists(EVENT_DICT_PATH)) {
        event_choices <- open_dataset(EVENT_DICT_PATH) %>%
          distinct(.data$outcome_name) %>%
          arrange(.data$outcome_name) %>%
          collect() %>%
          pull(.data$outcome_name)
        updateSelectizeInput(session, "event", choices = event_choices,
                             selected = "Myocardial infarction", server = TRUE)
      }
    })

    # Filter signals to selected pair
    selected_ts <- reactive({
      ds <- signals()
      req(ds)
      req(input$drug, input$event)
      ds %>%
        filter(.data$rxnorm_name == input$drug,
               .data$outcome_name == input$event) %>%
        collect() %>%
        arrange(.data$quarter)
    })

    output$timeline <- renderPlot({
      ts <- selected_ts()
      if (is.null(ts) || nrow(ts) == 0) {
        plot.new()
        graphics::text(0.5, 0.5, "No signals found for this pair.", cex = 1.5)
        return()
      }

      # Build caterpillar plot with base graphics + polygon for the ribbon
      ts <- mutate(ts,
        quarter_num = seq_along(.data$quarter),
        signal_tier = case_when(
          .data$eb05 >= 2 ~ "signal",
          .data$eb05 >= 1.5 ~ "watch",
          TRUE ~ "null"
        )
      )

      graphics::par(mar = c(4, 5, 3, 2))
      y_max <- max(c(ts$eb95, 3), na.rm = TRUE)
      y_min <- min(c(ts$eb05, 0.5), na.rm = TRUE)
      plot(
        ts$quarter_num, ts$eb50,
        type = "n",
        xlim = c(1, nrow(ts)),
        ylim = c(y_min, y_max),
        xaxt = "n",
        xlab = "Quarter",
        ylab = "Bayesian RR (EBGM with 95% CI)",
        main = paste0(input$drug, " + ", input$event),
        log = "y"
      )
      # Axis: quarter labels (every 4th or so to avoid clutter)
      step <- max(1, floor(nrow(ts) / 10))
      idx <- seq(1, nrow(ts), step)
      graphics::axis(1, at = idx, labels = ts$quarter[idx], las = 2, cex.axis = 0.8)

      # Reference lines
      graphics::abline(h = 1, col = "gray60", lty = 2)
      graphics::abline(h = 2, col = "firebrick", lty = 3)

      # Error bars (EB05 — EB95)
      graphics::segments(
        x0 = ts$quarter_num, y0 = ts$eb05,
        x1 = ts$quarter_num, y1 = ts$eb95,
        col = "steelblue", lwd = 2
      )
      # Point estimate
      color <- ifelse(ts$signal_tier == "signal", "firebrick",
                     ifelse(ts$signal_tier == "watch", "orange2", "steelblue"))
      graphics::points(ts$quarter_num, ts$eb50, pch = 19, col = color)
      # EWMA smoothed (if present)
      if ("ewma_eb05" %in% names(ts)) {
        graphics::lines(ts$quarter_num, ts$ewma_eb05, col = "darkgreen", lty = 1, lwd = 2)
      }

      graphics::legend(
        "topleft", bty = "n", cex = 0.9,
        legend = c("95% CI", "Point (EB50)", "Smoothed EB05", "Null (1)", "Signal (2)"),
        col = c("steelblue", "firebrick", "darkgreen", "gray60", "firebrick"),
        lty = c(1, NA, 1, 2, 3), lwd = c(2, NA, 2, 1, 1),
        pch = c(NA, 19, NA, NA, NA)
      )
    })
  })
}
