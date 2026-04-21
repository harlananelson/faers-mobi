# Signal Timeline — searchable / sortable DataTable of all signal-bearing
# (drug, event) pairs. Clicking a row renders a caterpillar plot of the
# credible intervals over time for that pair, plus a KNOWN/NOVEL badge
# against the drug's current FDA label.
#
# Reads precomputed signals from data/signals.parquet (produced offline by
# signal-compute). Never calls safetysignal live — all signals are
# precomputed on the local GPU box.

box::use(
  shiny[NS, moduleServer, tagList, plotOutput, renderPlot,
        req, reactive, reactiveVal, tags, div, h4, p, hr, fluidRow, column,
        uiOutput, renderUI, span, observeEvent, isolate],
  arrow[open_dataset, read_parquet],
  dplyr[filter, collect, pull, arrange, `%>%`, mutate, case_when,
        group_by, summarise, desc],
  DT[datatable, dataTableOutput, renderDataTable, formatRound, formatStyle, styleEqual],
)

SIGNALS_PATH <- "data/signals.parquet"
LABELS_PATH <- "data/fda_labels.parquet"
DIANA_PATH <- "data/diana_dictionary.parquet"

# ---- Novelty filter support ----
# Events that are medication errors, product-quality issues, or administrative
# reporting artefacts — not drug pharmacology. Dropped from "novel" candidates
# because calling them novel safety signals misrepresents the signal.
EVENT_BLACKLIST_EXACT <- c(
  "Intercepted drug dispensing error",
  "Drug dispensing error",
  "Medication error",
  "Wrong drug administered",
  "Product substitution issue",
  "Incorrect route of product administration",
  "Product dose omission issue",
  "Product quality issue",
  "Product solubility abnormal",
  "Product packaging issue",
  "Contamination product physical issue",
  "Expired product administered",
  "Accidental exposure to product",
  "Accidental overdose",
  "Off label use",
  "Drug abuse",
  "Drug dependence"
)
EVENT_BLACKLIST_PATTERNS <- c(
  "dispensing error", "administration error", "product.*error",
  "product use issue", "drug use error", "medication interaction"
)

.EVENT_STOP_WORDS <- c(
  "of", "the", "and", "in", "to", "a", "an", "on", "at", "with", "by",
  "for", "as", "is", "be", "from", "or", "under"
)

.event_in_label <- function(event, label_text) {
  ev <- tolower(event)
  lbl <- tolower(label_text)
  if (identical(lbl, "")) return(FALSE)
  if (grepl(ev, lbl, fixed = TRUE)) return(TRUE)
  words <- unlist(strsplit(ev, "[[:space:][:punct:]]+"))
  words <- words[nchar(words) >= 3 & !(words %in% .EVENT_STOP_WORDS)]
  if (length(words) == 0) return(FALSE)
  matched <- sum(vapply(words, function(w) grepl(w, lbl, fixed = TRUE), logical(1)))
  matched / length(words) >= 0.7
}

.event_is_blacklisted <- function(event) {
  if (event %in% EVENT_BLACKLIST_EXACT) return(TRUE)
  ev <- tolower(event)
  any(vapply(EVENT_BLACKLIST_PATTERNS,
             function(p) grepl(p, ev, ignore.case = TRUE), logical(1)))
}

# Resolve a raw drug name to an active-substance name using the DiAna
# dictionary (fusarolimichele/DiAna_package; ~348k raw-name -> substance
# pairs). Returns NA if unknown. Lookup is O(1) via match().
.diana_substance <- function(drug, diana) {
  if (is.null(diana) || nrow(diana) == 0) return(NA_character_)
  diana$substance[match(tolower(drug), diana$drugname)][1]
}

# Find the label row(s) for a drug name. Tries in order:
#   1. direct match on generic_name or brand_name (covers drugs without
#      a DiAna entry and rare/newly-approved products)
#   2. substance match: resolve the drug to its DiAna substance, then
#      filter labels by their precomputed substance column
# `lbl` is expected to carry a `substance` column (populated by the
# labels() reactive via DiAna when the dictionary is loaded).
.find_label_row <- function(lbl, drug, diana = NULL) {
  drug_lc <- tolower(drug)
  direct <- lbl[tolower(lbl$generic_name) == drug_lc |
                  (!is.na(lbl$brand_name) & tolower(lbl$brand_name) == drug_lc),
                , drop = FALSE]
  if (nrow(direct) > 0) return(direct)
  substance <- .diana_substance(drug, diana)
  if (is.na(substance) || !"substance" %in% names(lbl)) return(direct)
  lbl[!is.na(lbl$substance) & lbl$substance == substance, , drop = FALSE]
}

#' @export
ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        h4("Signals by drug and event"),
        p(
          "The top 2000 (drug, event) pairs by peak EB05, among those flagged ",
          "by \u22652 of 4 disproportionality methods (GPS/EBGM, PRR, ROR, IC). ",
          "The table is ",
          "searchable, sortable, and paginated \u2014 default sort is Peak EB05 ",
          "descending, so the strongest Bayesian signals are at the top. ",
          "Click any row to see the time-course plot and the label cross-check."
        ),
        p(tags$strong("Novel column:"),
          " \u201Cnovel\u201D = event (or \u226570% of its words) is absent from the ",
          "drug\u2019s boxed warning, contraindications, warnings, adverse ",
          "reactions, and (when cached) indications sections. Medication-error ",
          "and product-quality PTs are hidden. This is still a weak novelty ",
          "proxy \u2014 class effects and MedDRA-hierarchy synonyms are NOT ",
          "excluded. Treat novel rows as hypotheses to investigate, not ",
          "confirmed novel associations."),
        hr()
      )
    ),
    fluidRow(
      column(12, dataTableOutput(ns("signal_table")))
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

    # DiAna drug-name dictionary (raw-name -> active substance).
    # One-time load per session. ~348k rows, ~7 MB parquet.
    diana_dict <- reactive({
      if (!file.exists(DIANA_PATH)) return(NULL)
      read_parquet(DIANA_PATH)
    })

    # FDA label cache, augmented with a `substance` column derived from
    # DiAna: tries generic_name first, then brand_name. Lets both the
    # pair_stats novelty check and the KNOWN/NOVEL badge match signal
    # drugs to labels via DiAna-canonical substance, not raw text.
    labels <- reactive({
      if (!file.exists(LABELS_PATH)) return(NULL)
      lbl <- read_parquet(LABELS_PATH)
      diana <- diana_dict()
      if (!is.null(diana) && nrow(diana) > 0) {
        gn_sub <- diana$substance[match(tolower(lbl$generic_name), diana$drugname)]
        bn_sub <- diana$substance[match(tolower(lbl$brand_name), diana$drugname)]
        lbl$substance <- ifelse(!is.na(gn_sub), gn_sub, bn_sub)
      } else {
        lbl$substance <- NA_character_
      }
      lbl
    })

    # One-shot aggregation of every (drug, event) pair flagged by >=2 of 4
    # methods, with peak EWMA-smoothed EB05 and novelty flag. Used by the
    # datatable. Computed once per session (reactives cache).
    pair_stats <- reactive({
      ds <- signals()
      lbl <- labels()
      req(ds)
      # Aggregate on the arrow side — cheap, no full materialisation
      # Cap to the top 2000 by peak EB05. There are >140k pairs flagged by
      # >=2 methods; running the per-row label string match over all of
      # them blocks the page for minutes. The top 2000 covers the entire
      # clinically interesting range (weakest kept peak_eb05 will be far
      # below the signal threshold of 2).
      # Rows in signals where n_methods_flagged >= 2 are all signal-positive
      # by construction (is_signal_any requires >= 1 method). So the minimum
      # quarter over the filtered rows is the first quarter this pair was
      # flagged at the 2+ criterion — a useful "first signal" date.
      ps <- ds %>%
        filter(.data$n_methods_flagged >= 2) %>%
        group_by(.data$rxnorm_name, .data$outcome_name) %>%
        summarise(
          peak_eb05 = max(.data$ewma_eb05, na.rm = TRUE),
          n_methods_max = max(.data$n_methods_flagged, na.rm = TRUE),
          quarters_flagged = sum(.data$is_signal_any, na.rm = TRUE),
          first_signal = min(.data$quarter, na.rm = TRUE),
          latest_signal = max(.data$quarter, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(.data$peak_eb05)) %>%
        utils::head(2000) %>%
        collect()
      # Drop medication-error / admin / product-quality PTs — not drug effects
      ps <- ps[!vapply(ps$outcome_name, .event_is_blacklisted, logical(1)), , drop = FALSE]
      # Novelty column: TRUE (novel), FALSE (known), NA (no cached label)
      if (is.null(lbl)) {
        ps$novel <- NA
      } else {
        has_indications <- "indications_and_usage" %in% names(lbl)
        diana <- diana_dict()
        ps$novel <- mapply(function(drug, event) {
          row <- .find_label_row(lbl, drug, diana)
          if (nrow(row) == 0 || is.na(row$set_id[1])) return(NA)
          sections <- c(
            row$boxed_warning[1], row$contraindications[1],
            row$warnings_and_cautions[1], row$warnings[1], row$adverse_reactions[1]
          )
          if (has_indications) sections <- c(sections, row$indications_and_usage[1])
          sections[is.na(sections)] <- ""
          combined <- paste(sections, collapse = " \n ")
          !.event_in_label(event, combined)
        }, ps$rxnorm_name, ps$outcome_name)
      }
      ps
    })

    # Pick a reasonable default row so the plot has something to render on
    # first load. Falls back to row 1 if the preferred pair isn't present.
    .default_row <- function(ps) {
      preferred <- which(ps$rxnorm_name == "semaglutide" &
                         ps$outcome_name == "Cholecystitis acute")
      if (length(preferred) == 1) return(preferred)
      1L
    }

    output$signal_table <- renderDataTable({
      ps <- pair_stats()
      req(ps)
      display <- data.frame(
        Drug = ps$rxnorm_name,
        Event = ps$outcome_name,
        `Peak EB05` = round(ps$peak_eb05, 2),
        Methods = as.integer(ps$n_methods_max),
        Quarters = as.integer(ps$quarters_flagged),
        `First signal` = ps$first_signal,
        `Latest signal` = ps$latest_signal,
        Novel = ifelse(is.na(ps$novel), "?", ifelse(ps$novel, "novel", "known")),
        check.names = FALSE
      )
      datatable(
        display,
        selection = list(mode = "single", selected = .default_row(ps)),
        rownames = FALSE,
        filter = "top",
        options = list(
          pageLength = 25,
          lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
          order = list(list(2, "desc")),
          searchHighlight = TRUE,
          columnDefs = list(list(className = "dt-right", targets = c(2, 3, 4)))
        )
      ) |>
        formatStyle(
          "Novel",
          backgroundColor = styleEqual(
            c("novel", "known", "?"),
            c("#ffe8e8", "#e8f5e8", "#f0f0f0")
          )
        )
    })

    # Which (drug, event) the user has selected from the table. Starts from
    # the default row pre-selected by renderDataTable above.
    selected_pair <- reactive({
      ps <- pair_stats()
      req(ps)
      i <- input$signal_table_rows_selected
      if (is.null(i) || length(i) == 0) i <- .default_row(ps)
      list(drug = ps$rxnorm_name[i], event = ps$outcome_name[i])
    })

    # Filter signals to selected pair for plotting
    selected_ts <- reactive({
      ds <- signals()
      sp <- selected_pair()
      req(ds, sp$drug, sp$event)
      ds %>%
        filter(.data$rxnorm_name == sp$drug,
               .data$outcome_name == sp$event) %>%
        collect() %>%
        arrange(.data$quarter)
    })

    # Known-vs-novel badge: does the event appear in the drug's current label?
    output$known_badge <- renderUI({
      sp <- selected_pair()
      req(sp$drug, sp$event)
      lbl_df <- labels()
      if (is.null(lbl_df)) {
        return(tags$div(class = "alert alert-secondary small py-2",
                        "Label cross-reference not loaded."))
      }
      row <- .find_label_row(lbl_df, sp$drug, diana_dict())
      if (nrow(row) == 0 || is.na(row$set_id[1])) {
        return(tags$div(class = "alert alert-secondary small py-2",
                        tags$strong("No label cached"),
                        " for ", tags$em(sp$drug), " \u2014 can't determine if this signal is known."))
      }
      nz <- function(x) if (is.null(x) || is.na(x)) "" else x
      sections <- c(
        Boxed = nz(row$boxed_warning[1]),
        Contraindications = nz(row$contraindications[1]),
        `Warnings/Precautions` = paste(nz(row$warnings_and_cautions[1]), nz(row$warnings[1])),
        `Adverse Reactions` = nz(row$adverse_reactions[1])
      )
      hits <- names(sections)[vapply(sections, function(s) .event_in_label(sp$event, s), logical(1))]
      if (length(hits) == 0) {
        tags$div(class = "alert alert-danger small py-2",
                 tags$strong("NOVEL"),
                 " \u2014 \"", tags$em(sp$event), "\" is not mentioned in ",
                 tags$em(sp$drug), "'s current FDA label. Treat as a hypothesis worth investigating.")
      } else {
        priority <- c("Boxed", "Contraindications", "Warnings/Precautions", "Adverse Reactions")
        top <- priority[priority %in% hits][1]
        tags$div(class = "alert alert-success small py-2",
                 tags$strong("KNOWN"),
                 " \u2014 \"", tags$em(sp$event), "\" appears in ",
                 tags$em(sp$drug), "'s label (",
                 paste(hits, collapse = ", "), "). This signal is already documented.")
      }
    })

    output$timeline <- renderPlot({
      ts <- selected_ts()
      sp <- selected_pair()
      if (is.null(ts) || nrow(ts) == 0) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No signals found for this pair.", cex = 1.5)
        return()
      }

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
        main = paste0(sp$drug, " + ", sp$event),
        log = "y"
      )
      step <- max(1, floor(nrow(ts) / 10))
      idx <- seq(1, nrow(ts), step)
      graphics::axis(1, at = idx, labels = ts$quarter[idx], las = 2, cex.axis = 0.8)

      graphics::abline(h = 1, col = "gray60", lty = 2)
      graphics::abline(h = 2, col = "firebrick", lty = 3)

      graphics::segments(
        x0 = ts$quarter_num, y0 = ts$eb05,
        x1 = ts$quarter_num, y1 = ts$eb95,
        col = "steelblue", lwd = 2
      )
      color <- ifelse(ts$signal_tier == "signal", "firebrick",
                     ifelse(ts$signal_tier == "watch", "orange2", "steelblue"))
      graphics::points(ts$quarter_num, ts$eb50, pch = 19, col = color)
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
