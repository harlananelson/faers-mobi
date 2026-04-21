# Signal Timeline — caterpillar plot of credible intervals over time for a
# single (drug, event) pair.
#
# Reads precomputed signals from /srv/shiny-server/faers-mobi/data/signals.parquet
# (produced offline by signal-compute). Never calls safetysignal live — all
# signals are precomputed on the local GPU box.

box::use(
  shiny[NS, moduleServer, tagList, selectizeInput, plotOutput, renderPlot,
        req, reactive, reactiveVal, tags, div, h4, p, hr, fluidRow, column,
        withProgress, updateSelectizeInput, uiOutput, renderUI, span,
        observeEvent, isolate, tableOutput, renderTable],
  arrow[open_dataset, read_parquet],
  dplyr[filter, collect, pull, arrange, distinct, `%>%`, mutate, case_when,
        group_by, summarise, desc],
)

SIGNALS_PATH <- "data/signals.parquet"
LABELS_PATH <- "data/fda_labels.parquet"

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
      column(12,
        h4("Top 20 novel signals — not in FDA label"),
        p("Computed across all pairs flagged by \u22653 of 4 methods, ranked by ",
          "peak EWMA-smoothed EB05. \u201CNovel\u201D means the event is not ",
          "mentioned in the drug\u2019s current FDA label (limited to drugs ",
          "with cached openFDA label data). These are hypotheses worth ",
          "investigating, not confirmed associations."),
        tableOutput(ns("top_novel_table")),
        hr()
      )
    ),
    fluidRow(
      column(6, selectizeInput(ns("drug"), "Drug (RxNorm):",
                                choices = NULL, multiple = FALSE)),
      column(6, selectizeInput(ns("event"), "Event (MedDRA PT):",
                                choices = NULL, multiple = FALSE))
    ),
    fluidRow(column(12, uiOutput(ns("known_badge")))),
    fluidRow(column(12, plotOutput(ns("timeline"), height = "500px"))),
    fluidRow(
      column(12,
        tags$small(tags$em(
          "Disclaimer: disproportionate reporting is a statistical pattern, ",
          "not evidence of causation. Signals are hypotheses requiring further ",
          "investigation. 'Known' means the event appears in the drug's current ",
          "FDA label; 'Novel' means it does not (label coverage is limited to ",
          "drugs we have cached openFDA label data for)."
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

    # One-shot cache of distinct (drug, event) pairs that actually have signal
    # rows. Drives both the initial dropdown population and the cross-filter
    # observers below — users can only pick combinations where data exists.
    pairs_cache <- reactive({
      ds <- signals()
      if (is.null(ds)) return(NULL)
      ds %>%
        distinct(.data$rxnorm_name, .data$outcome_name) %>%
        collect()
    })

    # Initial population: both dropdowns from pairs_cache, defaulting to the
    # semaglutide / Cholecystitis acute pair (clear GLP-1 signal) when available.
    shiny::observe({
      pairs <- pairs_cache()
      req(pairs)
      drugs <- sort(unique(pairs$rxnorm_name))
      events <- sort(unique(pairs$outcome_name))
      default_drug <- if ("semaglutide" %in% drugs) "semaglutide" else drugs[1]
      default_event <- if ("Cholecystitis acute" %in% events) "Cholecystitis acute" else events[1]
      updateSelectizeInput(session, "drug", choices = drugs,
                           selected = default_drug, server = TRUE)
      updateSelectizeInput(session, "event", choices = events,
                           selected = default_event, server = TRUE)
    })

    # Flags to break the cross-filter loop. When an observer programmatically
    # updates the other select, it sets the partner's skip flag. The partner
    # observer fires (because updateSelectizeInput re-emits the input),
    # consumes the flag, and returns without running the filter.
    skip_next_drug <- reactiveVal(FALSE)
    skip_next_event <- reactiveVal(FALSE)

    # When the drug changes, narrow the event list to events that have a
    # signal row for the selected drug. Clearing the drug restores the full
    # list. Keeps the current event selected if still valid.
    observeEvent(input$drug, {
      if (isolate(skip_next_drug())) {
        skip_next_drug(FALSE)
        return()
      }
      pairs <- pairs_cache()
      req(pairs)
      events <- if (is.null(input$drug) || !nzchar(input$drug)) {
        sort(unique(pairs$outcome_name))
      } else {
        pairs %>% filter(.data$rxnorm_name == input$drug) %>%
          pull(.data$outcome_name) %>% unique() %>% sort()
      }
      current_event <- isolate(input$event)
      keep <- if (!is.null(current_event) && nzchar(current_event) && current_event %in% events) current_event else character(0)
      skip_next_event(TRUE)
      updateSelectizeInput(session, "event", choices = events, selected = keep, server = TRUE)
    }, ignoreInit = TRUE)

    # Symmetric: when the event changes, narrow the drug list.
    observeEvent(input$event, {
      if (isolate(skip_next_event())) {
        skip_next_event(FALSE)
        return()
      }
      pairs <- pairs_cache()
      req(pairs)
      drugs <- if (is.null(input$event) || !nzchar(input$event)) {
        sort(unique(pairs$rxnorm_name))
      } else {
        pairs %>% filter(.data$outcome_name == input$event) %>%
          pull(.data$rxnorm_name) %>% unique() %>% sort()
      }
      current_drug <- isolate(input$drug)
      keep <- if (!is.null(current_drug) && nzchar(current_drug) && current_drug %in% drugs) current_drug else character(0)
      skip_next_drug(TRUE)
      updateSelectizeInput(session, "drug", choices = drugs, selected = keep, server = TRUE)
    }, ignoreInit = TRUE)

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

    # Top 20 pairs flagged by >=3 of 4 methods whose event is NOT mentioned
    # in the drug's current FDA label. Computed once per session.
    top_novel <- reactive({
      ds <- signals()
      lbl <- labels()
      req(ds, lbl)
      top <- ds %>%
        filter(.data$n_methods_flagged >= 3) %>%
        group_by(.data$rxnorm_name, .data$outcome_name) %>%
        summarise(
          peak_eb05 = max(.data$ewma_eb05, na.rm = TRUE),
          n_methods_max = max(.data$n_methods_flagged, na.rm = TRUE),
          quarters_flagged = sum(.data$is_signal_any, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(.data$peak_eb05)) %>%
        utils::head(200) %>%
        collect()
      top$is_novel <- mapply(function(drug, event) {
        row <- lbl[lbl$generic_name == tolower(drug), , drop = FALSE]
        if (nrow(row) == 0 || is.na(row$set_id[1])) return(NA)
        ev <- tolower(event)
        sections <- c(
          row$boxed_warning[1], row$contraindications[1],
          row$warnings_and_cautions[1], row$warnings[1], row$adverse_reactions[1]
        )
        sections[is.na(sections)] <- ""
        !any(grepl(ev, tolower(sections), fixed = TRUE))
      }, top$rxnorm_name, top$outcome_name)
      top <- top[!is.na(top$is_novel) & top$is_novel, , drop = FALSE]
      top$is_novel <- NULL
      utils::head(top, 20)
    })

    output$top_novel_table <- renderTable({
      tn <- top_novel()
      req(tn)
      if (nrow(tn) == 0) {
        return(data.frame(Message = "No novel signals at current criteria."))
      }
      out <- data.frame(
        Drug = tn$rxnorm_name,
        Event = tn$outcome_name,
        `Peak EB05` = round(tn$peak_eb05, 2),
        `Methods (max of 4)` = as.integer(tn$n_methods_max),
        `Quarters flagged` = as.integer(tn$quarters_flagged),
        check.names = FALSE
      )
      out
    }, digits = 2, striped = TRUE, hover = TRUE, spacing = "s")

    # FDA label cache (one-time load per session)
    labels <- reactive({
      if (!file.exists(LABELS_PATH)) return(NULL)
      read_parquet(LABELS_PATH)
    })

    # Known-vs-novel badge: does the event appear in the drug's current label?
    output$known_badge <- renderUI({
      req(input$drug, input$event)
      lbl_df <- labels()
      if (is.null(lbl_df)) {
        return(tags$div(class = "alert alert-secondary small py-2",
                        "Label cross-reference not loaded."))
      }
      row <- lbl_df[lbl_df$generic_name == tolower(input$drug), , drop = FALSE]
      if (nrow(row) == 0 || is.na(row$set_id[1])) {
        return(tags$div(class = "alert alert-secondary small py-2",
                        tags$strong("No label cached"),
                        " for ", tags$em(input$drug), " — can't determine if this signal is known."))
      }
      ev <- tolower(input$event)
      # Check each section; report the most serious hit
      nz <- function(x) if (is.null(x) || is.na(x)) "" else x
      sections <- c(
        Boxed = nz(row$boxed_warning[1]),
        Contraindications = nz(row$contraindications[1]),
        `Warnings/Precautions` = paste(nz(row$warnings_and_cautions[1]), nz(row$warnings[1])),
        `Adverse Reactions` = nz(row$adverse_reactions[1])
      )
      hits <- names(sections)[vapply(sections, function(s) grepl(ev, tolower(s), fixed = TRUE), logical(1))]
      if (length(hits) == 0) {
        tags$div(class = "alert alert-danger small py-2",
                 tags$strong("NOVEL"),
                 " — \"", tags$em(input$event), "\" is not mentioned in ",
                 tags$em(input$drug), "'s current FDA label. Treat as a hypothesis worth investigating.")
      } else {
        priority <- c("Boxed", "Contraindications", "Warnings/Precautions", "Adverse Reactions")
        top <- priority[priority %in% hits][1]
        tags$div(class = "alert alert-success small py-2",
                 tags$strong("KNOWN"),
                 " — \"", tags$em(input$event), "\" appears in ",
                 tags$em(input$drug), "'s label (",
                 paste(hits, collapse = ", "), "). This signal is already documented.")
      }
    })

    output$timeline <- renderPlot({
      ts <- selected_ts()
      if (is.null(ts) || nrow(ts) == 0) {
        graphics::plot.new()
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
      graphics::plot(
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
