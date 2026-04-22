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
FIRST_APPROVAL_PATH <- "data/first_approval.parquet"
MEDDRA_PATH <- "data/meddra_hierarchy.parquet"
ATC_PATH <- "data/atc_classes.parquet"
BLACKLIST_EXACT_PATH <- "data/event_blacklist_exact.csv"
BLACKLIST_PATTERN_PATH <- "data/event_blacklist_patterns.csv"
TRIAGE_PATH <- "data/signal_triage.csv"
WATCHLIST_PATH <- "data/signal_watchlist.csv"

# Load a control file's single column. Returns character(0) when the
# file is missing or the column doesn't exist, so the rest of the
# pipeline is unaffected by a broken/missing control file.
.load_control_column <- function(path, col) {
  if (!file.exists(path)) return(character(0))
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8"),
    error = function(e) NULL
  )
  if (is.null(df) || !col %in% names(df)) return(character(0))
  vals <- as.character(df[[col]])
  vals[nchar(vals) > 0]
}

# ---- Novelty filter support ----
# Events that are medication errors, product-quality issues, or administrative
# reporting artefacts — not drug pharmacology. Dropped from "novel" candidates
# because calling them novel safety signals misrepresents the signal.
# The control files at data/event_blacklist_*.csv are authoritative; the
# in-file defaults below are a fallback when those files are missing.
.FALLBACK_BLACKLIST_EXACT <- c(
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
  "Drug dependence",
  # Added 2026-04-21 from top-20 novel audit (AI/reviews/top20-novel-verification.md).
  # These are all MedDRA SOC "Product issues" or medication-error PTs dressed
  # as AEs — not drug pharmacology.
  "Product contamination",
  "Product contamination physical",
  "Product contamination microbial",
  "Product contamination chemical",
  "Product cleaning inadequate",
  "Product residue present",
  "Product deposit",
  "Labelled drug-disease interaction medication error",
  "Product physical issue",
  "Injection site extravasation",
  # Added 2026-04-22 from Grok top-100 blacklist review
  # (AI/plans/blacklist.md). Clear-noise PTs only.
  "Out of specification product use",
  "Medication residue present",
  "Counterfeit product administered",
  "Reaction to excipient",
  "Food interaction",
  "Papillary deformity",
  "Foreign body",
  "Foreign body in respiratory tract",
  "Removal of foreign body",
  "Contrast media reaction",
  "Contrast media toxicity",
  "Fibrin",
  "Shunt stenosis",
  "Shunt occlusion",
  "Retinal exudates",
  "Enterobacter test positive",
  "Burkholderia test positive",
  "Burkholderia infection",
  "Culture positive",
  "Sputum culture positive",
  "Infantile spitting up",
  "Pathogen resistance",
  "Expulsion of medication"
)
.FALLBACK_BLACKLIST_PATTERNS <- c(
  "dispensing error", "administration error", "product.*error",
  "product use issue", "drug use error", "medication interaction",
  # Broaden: any PT starting with "product " is SOC Product issues
  "^product ",
  # Any PT ending in " medication error" or " dispensing error"
  "medication error$",
  # Drug-disease interaction errors
  "drug-disease interaction",
  # Drug preparation / storage / administration issues
  "drug preparation error",
  "product storage error",
  # 2026-04-22 additions from Grok top-100 review
  "test positive$",
  "antigen increased$",
  "^foreign body",
  "^contrast media"
)

# Load the canonical blacklists from the control files. Falls back to the
# hardcoded defaults if the CSVs aren't present on disk.
.loaded_exact <- .load_control_column(BLACKLIST_EXACT_PATH, "pt")
.loaded_patt  <- .load_control_column(BLACKLIST_PATTERN_PATH, "pattern")
EVENT_BLACKLIST_EXACT <- if (length(.loaded_exact) > 0) .loaded_exact else .FALLBACK_BLACKLIST_EXACT
EVENT_BLACKLIST_PATTERNS <- if (length(.loaded_patt) > 0) .loaded_patt else .FALLBACK_BLACKLIST_PATTERNS

.EVENT_STOP_WORDS <- c(
  "of", "the", "and", "in", "to", "a", "an", "on", "at", "with", "by",
  "for", "as", "is", "be", "from", "or", "under"
)

# Normalize British <-> American medical spellings so that a US label
# that says "anemia / leukemia / tumor" matches a MedDRA PT written
# in British spelling ("anaemia / leukaemia / tumour") and vice-versa.
.normalize_spelling <- function(x) {
  # British -> American medical spellings (one-way, American canonical).
  x <- gsub("aemia",   "emia",   x, fixed = TRUE)
  x <- gsub("oedema",  "edema",  x, fixed = TRUE)
  x <- gsub("oesoph",  "esoph",  x, fixed = TRUE)
  x <- gsub("haema",   "hema",   x, fixed = TRUE)
  x <- gsub("haemo",   "hemo",   x, fixed = TRUE)
  x <- gsub("tumour",  "tumor",  x, fixed = TRUE)
  x <- gsub("leukaem", "leukem", x, fixed = TRUE)
  x <- gsub("paediat", "pediat", x, fixed = TRUE)
  x <- gsub("diarrhoea", "diarrhea", x, fixed = TRUE)
  x <- gsub("colour",  "color",  x, fixed = TRUE)  # discolouration / colour
  x <- gsub("fibres",  "fibers", x, fixed = TRUE)
  x <- gsub("fibre",   "fiber",  x, fixed = TRUE)
  x <- gsub("centre",  "center", x, fixed = TRUE)
  x <- gsub("metre",   "meter",  x, fixed = TRUE)
  # Curated clinical synonyms where UMLS atoms often miss the pairing.
  # Each entry picks a canonical form and rewrites to it on both sides.
  x <- gsub("adrenocortical", "adrenal",     x, fixed = TRUE)
  x <- gsub("lymphoblastic",  "lymphocytic", x, fixed = TRUE)  # MedDRA canonical
  x <- gsub("relapsed",       "recurrent",   x, fixed = TRUE)
  x <- gsub("refractory",     "recurrent",   x, fixed = TRUE)
  x <- gsub("staining",       "discoloration", x, fixed = TRUE) # tooth/skin staining
  x
}

# Fuzzy match of an event string against arbitrary label text.
# Case-insensitive; strips stop-words and short tokens; matches when
# either (a) the event is a direct substring of the label or (b) at
# least `threshold` fraction of event's non-trivial words appear
# somewhere in the label. Spelling normalized on both sides before
# comparison. Threshold default 0.7 for direct event matches;
# callers can pass 0.6 for MedDRA-synonym matches where some semantic
# slack is appropriate.
.event_in_label <- function(event, label_text, threshold = 0.7) {
  ev <- .normalize_spelling(tolower(event))
  lbl <- .normalize_spelling(tolower(label_text))
  if (identical(lbl, "")) return(FALSE)
  if (grepl(ev, lbl, fixed = TRUE)) return(TRUE)
  words <- unlist(strsplit(ev, "[[:space:][:punct:]]+"))
  words <- words[nchar(words) >= 3 & !(words %in% .EVENT_STOP_WORDS)]
  if (length(words) == 0) return(FALSE)
  matched <- sum(vapply(words, function(w) grepl(w, lbl, fixed = TRUE), logical(1)))
  matched / length(words) >= threshold
}

.event_is_blacklisted <- function(event) {
  if (event %in% EVENT_BLACKLIST_EXACT) return(TRUE)
  ev <- tolower(event)
  any(vapply(EVENT_BLACKLIST_PATTERNS,
             function(p) grepl(p, ev, ignore.case = TRUE), logical(1)))
}

# Expanded in-label check: returns TRUE if the event, any of its
# non-trivial word-stems, OR any MedDRA hierarchy synonym (via UMLS CUI)
# is present in the label text. Covers cases where the label describes
# the same concept with non-MedDRA wording (e.g. PT "Spinal cord
# infarction" vs label "spinal cord ischaemia").
.event_in_label_expanded <- function(event, label_text, meddra_row = NULL) {
  if (.event_in_label(event, label_text, threshold = 0.7)) return(TRUE)
  if (is.null(meddra_row) || nrow(meddra_row) == 0) return(FALSE)
  syns <- meddra_row$syns_list[[1]]
  syns <- syns[nchar(syns) >= 5]
  if (length(syns) == 0) return(FALSE)
  # Each synonym gets a slightly looser word-level match (0.6) against
  # the label. Synonyms are already semantically equivalent to the
  # event, so giving them 10 percentage points more slack is sound.
  any(vapply(syns, function(s) .event_in_label(s, label_text, threshold = 0.6),
             logical(1)))
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
        tags$div(class = "alert alert-info small py-2",
          tags$strong("Note:"),
          " the first load takes 10-20 seconds while the app computes label ",
          "matches, MedDRA synonyms, and class stats for the top 2000 pairs. ",
          "Subsequent interactions are fast."),
        p(
          "The top 2000 (drug, event) pairs by peak EB05, among those flagged ",
          "by \u22652 of 4 disproportionality methods (GPS/EBGM, PRR, ROR, IC). ",
          "The table is searchable, sortable, and paginated. Default view: ",
          tags$strong("Novel"), " pairs with \u22653 quarters of signal, sorted by ",
          tags$strong("Adj EB05"), " descending (peak EB05 after Weber-effect ",
          "shrinkage for drugs <5 years on market). Click any row to see the ",
          "time-course plot and the label cross-check."
        ),
        p(tags$strong("Novel column:"),
          " \u201Cnovel\u201D means the event is absent from the drug\u2019s boxed ",
          "warning, contraindications, warnings, adverse reactions, and ",
          "indications sections \u2014 after MedDRA-synonym expansion (UMLS CUI) ",
          "and British\u2194American + clinical-term normalization (anaemia/anemia, ",
          "adrenocortical/adrenal, lymphoblastic/lymphocytic, etc). Medication-",
          "error, product-quality, and administration PTs are hidden. ",
          tags$strong("Class co-flags"),
          " = number of other drugs in the same ATC4 class that also flag this ",
          "event in the top-2000 slice (1 = drug-specific; >1 suggests class ",
          "effect). These filters substantially reduce false-positive ",
          "\u201Cnovel\u201D flags, but treat remaining rows as hypotheses to ",
          "investigate, not confirmed novel associations."),
        hr()
      )
    ),
    fluidRow(
      column(12, dataTableOutput(ns("signal_table")))
    ),
    fluidRow(column(12, uiOutput(ns("known_badge")))),
    fluidRow(column(12, uiOutput(ns("event_description")))),
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

    # FDA Orange Book first-approval-for-marketing dates, keyed by active
    # substance. ~2100 single-ingredient approved drugs.
    first_approval <- reactive({
      if (!file.exists(FIRST_APPROVAL_PATH)) return(NULL)
      read_parquet(FIRST_APPROVAL_PATH)
    })

    # MedDRA PT hierarchy via UMLS: for each top-signal PT we cache its
    # CUI and all same-CUI synonyms across vocabularies. Used by the
    # Novel check so that label text using non-MedDRA wording still
    # resolves to "known".
    meddra <- reactive({
      if (!file.exists(MEDDRA_PATH)) return(NULL)
      df <- read_parquet(MEDDRA_PATH)
      df$syns_list <- strsplit(df$synonyms, ";", fixed = TRUE)
      df
    })

    # ATC class binder from DiAna OSF: substance -> atc_code +
    # atc_class4 (chemical subgroup, e.g. "Proton pump inhibitors").
    # Used for the Class column in the DT; enables class-effect
    # visual scanning.
    atc_classes <- reactive({
      if (!file.exists(ATC_PATH)) return(NULL)
      read_parquet(ATC_PATH)
    })

    # Manual-triage control file: curated classifications for
    # specific (drug, event) pairs with explanations. Edit the CSV to
    # add / revise entries; changes take effect on next app restart.
    # Schema: drug, event, classification, explanation, reviewer, date.
    triage <- reactive({
      if (!file.exists(TRIAGE_PATH)) return(NULL)
      df <- tryCatch(
        utils::read.csv(TRIAGE_PATH, stringsAsFactors = FALSE, encoding = "UTF-8"),
        error = function(e) NULL
      )
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$drug_lc <- tolower(trimws(df$drug))
      df$event_lc <- tolower(trimws(df$event))
      df
    })

    # Watchlist control file: pairs to track over time regardless of
    # triage status. Schema: drug, event, reason, added_by, added_date.
    watchlist <- reactive({
      if (!file.exists(WATCHLIST_PATH)) return(NULL)
      df <- tryCatch(
        utils::read.csv(WATCHLIST_PATH, stringsAsFactors = FALSE, encoding = "UTF-8"),
        error = function(e) NULL
      )
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$drug_lc <- tolower(trimws(df$drug))
      df$event_lc <- tolower(trimws(df$event))
      df
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
      # Attach DiAna-resolved substance and first-approval date per drug
      diana <- diana_dict()
      appr <- first_approval()
      if (!is.null(diana) && nrow(diana) > 0) {
        ps$substance <- diana$substance[match(tolower(ps$rxnorm_name), diana$drugname)]
      } else {
        ps$substance <- NA_character_
      }
      if (!is.null(appr) && nrow(appr) > 0) {
        ps$first_approval <- appr$first_approval[match(ps$substance, appr$substance)]
        need_fallback <- is.na(ps$first_approval)
        ps$first_approval[need_fallback] <-
          appr$first_approval[match(tolower(ps$rxnorm_name[need_fallback]),
                                    appr$substance)]
      } else {
        ps$first_approval <- as.Date(NA)
      }
      # Weber-effect adjustment: shrink peak_eb05 for recently marketed
      # drugs. New drugs with tiny background counts produce inflated EB05;
      # this post-hoc correction (eb05 * (1 - 0.6 * exp(-years_on_market)))
      # pulls the score down ~36% at 6 months, ~22% at 1 year, ~0% at 5+.
      # Drugs without a known approval date get no adjustment (assumed
      # established). See NOVELTY_FILTER_ROADMAP.md for the proper adaptive
      # prior path (pipeline-side).
      today <- Sys.Date()
      ps$years_on_market <- as.numeric(today - ps$first_approval) / 365.25
      shrink <- 1 - 0.6 * exp(-pmax(ps$years_on_market, 0))
      ps$adj_eb05 <- ifelse(is.na(shrink), ps$peak_eb05, ps$peak_eb05 * shrink)
      # ATC class4 (chemical subgroup) via DiAna substance
      atc <- atc_classes()
      if (!is.null(atc) && nrow(atc) > 0) {
        ps$atc_class <- atc$atc_class4[match(ps$substance, atc$substance)]
      } else {
        ps$atc_class <- NA_character_
      }
      # Class co-flags: for each (class, event), count how many distinct
      # drugs in the same ATC4 class also flag this event (among top-2000).
      # A high number means "class effect" (likely already labeled at the
      # class level even if missing from a specific drug's text);
      # 1 means drug-specific.
      has_class <- !is.na(ps$atc_class) & !is.na(ps$outcome_name)
      ps$class_co_flags <- 1L
      if (any(has_class)) {
        key <- paste(ps$atc_class, ps$outcome_name, sep = "\x1f")
        counts <- table(key[has_class])
        ps$class_co_flags[has_class] <-
          as.integer(counts[match(key[has_class], names(counts))])
      }
      # Join manual-triage classifications. Match on lowercased drug+event.
      tri <- triage()
      key <- paste(tolower(ps$rxnorm_name), tolower(ps$outcome_name),
                   sep = "\x1f")
      if (!is.null(tri)) {
        tri_key <- paste(tri$drug_lc, tri$event_lc, sep = "\x1f")
        idx <- match(key, tri_key)
        ps$triage <- ifelse(is.na(idx), "", tri$classification[idx])
      } else {
        ps$triage <- ""
      }
      # Join watchlist. Watched pairs are marked with a star.
      wl <- watchlist()
      if (!is.null(wl)) {
        wl_key <- paste(wl$drug_lc, wl$event_lc, sep = "\x1f")
        ps$watch <- ifelse(key %in% wl_key, "\u2605", "")
      } else {
        ps$watch <- ""
      }
      # Novelty column: TRUE (novel), FALSE (known), NA (no cached label)
      if (is.null(lbl)) {
        ps$novel <- NA
        ps$treats <- NA
      } else {
        has_indications <- "indications_and_usage" %in% names(lbl)
        mh <- meddra()
        results <- mapply(function(drug, event) {
          row <- .find_label_row(lbl, drug, diana)
          if (nrow(row) == 0 || is.na(row$set_id[1])) return(c(NA, NA))
          sections <- c(
            row$boxed_warning[1], row$contraindications[1],
            row$warnings_and_cautions[1], row$warnings[1], row$adverse_reactions[1]
          )
          if (has_indications) sections <- c(sections, row$indications_and_usage[1])
          sections[is.na(sections)] <- ""
          combined <- paste(sections, collapse = " \n ")
          mrow <- if (!is.null(mh)) mh[mh$pt == event, , drop = FALSE] else NULL
          novel <- !.event_in_label_expanded(event, combined, mrow)
          treats <- FALSE
          if (has_indications) {
            ind <- if (is.na(row$indications_and_usage[1])) "" else row$indications_and_usage[1]
            treats <- .event_in_label_expanded(event, ind, mrow)
          }
          c(novel, treats)
        }, ps$rxnorm_name, ps$outcome_name)
        ps$novel <- as.logical(results[1, ])
        ps$treats <- as.logical(results[2, ])
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
      # Drug label: show "brand (ingredient)" when DiAna resolved them to
      # different strings, else just the raw name. Lets users see both the
      # name they know and the canonical active ingredient.
      drug_col <- ifelse(
        !is.na(ps$substance) & tolower(ps$rxnorm_name) != ps$substance,
        paste0(ps$rxnorm_name, " (", ps$substance, ")"),
        ps$rxnorm_name
      )
      display <- data.frame(
        Watch = ps$watch,
        Drug = drug_col,
        Event = ps$outcome_name,
        `Peak EB05` = round(ps$peak_eb05, 2),
        `Adj EB05` = round(ps$adj_eb05, 2),
        Quarters = as.integer(ps$quarters_flagged),
        `First FDA Report` = ps$first_signal,
        `Latest Report` = ps$latest_signal,
        `Approval Year` = ifelse(is.na(ps$first_approval), "",
                                 format(ps$first_approval, "%Y")),
        `Yrs on Market` = ifelse(is.na(ps$first_approval), NA_integer_,
                                  as.integer(floor(ps$years_on_market))),
        Class = ifelse(is.na(ps$atc_class), "", ps$atc_class),
        `Class co-flags` = as.integer(ps$class_co_flags),
        Treats = ifelse(is.na(ps$treats), "", ifelse(ps$treats, "yes", "")),
        Triage = ps$triage,
        Novel = ifelse(is.na(ps$novel), "?", ifelse(ps$novel, "novel", "known")),
        check.names = FALSE
      )
      # Column indices (0-based): 0 Watch, 1 Drug, 2 Event, 3 Peak EB05,
      # 4 Adj EB05, 5 Quarters, 6 First FDA Report, 7 Latest Report,
      # 8 Approval Year, 9 Yrs on Market, 10 Class, 11 Class co-flags,
      # 12 Treats, 13 Triage, 14 Novel. Default sort: Adj EB05 desc
      # (col 4). Default filter: Novel = "novel" and Quarters >= 3.
      datatable(
        display,
        selection = list(mode = "single", selected = .default_row(ps)),
        rownames = FALSE,
        filter = "top",
        extensions = "Buttons",
        options = list(
          pageLength = 25,
          lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
          order = list(list(4, "desc")),
          searchHighlight = TRUE,
          dom = 'Blfrtip',
          buttons = list(list(extend = "csv", text = "Download CSV",
                              filename = "signals")),
          searchCols = list(
            NULL, NULL, NULL, NULL, NULL,
            list(search = "3 ... 9999"),
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
            list(search = "novel")
          ),
          columnDefs = list(list(className = "dt-right",
                                 targets = c(3, 4, 5, 9, 11)))
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

    output$event_description <- renderUI({
      sp <- selected_pair()
      req(sp$event)
      mh <- meddra()
      if (is.null(mh) || !"definition" %in% names(mh)) return(NULL)
      row <- mh[mh$pt == sp$event, , drop = FALSE]
      if (nrow(row) == 0 || is.na(row$definition[1]) || nchar(row$definition[1]) == 0) {
        return(NULL)
      }
      tags$div(
        class = "alert alert-secondary small py-2",
        tags$strong("About the event: "),
        tags$em(sp$event), " \u2014 ", row$definition[1]
      )
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
      # Connect EB50 point estimates with a straight line (no EWMA
      # smoothing — the points are already Bayes-shrunk posterior
      # medians, so smoothing them is redundant).
      graphics::lines(ts$quarter_num, ts$eb50, col = "gray40", lwd = 1, lty = 1)
      color <- ifelse(ts$signal_tier == "signal", "firebrick",
                     ifelse(ts$signal_tier == "watch", "orange2", "steelblue"))
      graphics::points(ts$quarter_num, ts$eb50, pch = 19, col = color)

      graphics::legend(
        "topleft", bty = "n", cex = 0.9,
        legend = c("95% CI", "Point (EB50)", "Null (1)", "Signal (2)"),
        col = c("steelblue", "firebrick", "gray60", "firebrick"),
        lty = c(1, NA, 2, 3), lwd = c(2, NA, 1, 1),
        pch = c(NA, 19, NA, NA)
      )
    })
  })
}
