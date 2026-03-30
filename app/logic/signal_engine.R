# Signal detection engine using safetysignal package

box::use(
  safetysignal[compute_observed_expected, fit_prior, compute_posterior,
               posterior_percentile, detect_signals],
  dplyr[filter, arrange, desc, select, mutate],
  tibble[tibble],
)

#' Run GPS signal detection on drug-event data
#'
#' @param data A tibble with columns: product, event, and optionally count.
#' @param threshold EB05 threshold for signal detection.
#' @return A tibble with signal detection results.
#' @export
run_signal_detection <- function(data, threshold = 1) {
  oe <- compute_observed_expected(data, product, event)

  if (nrow(oe) < 10) {
    return(oe)
  }

  prior <- fit_prior(oe)

  oe |>
    compute_posterior(prior) |>
    posterior_percentile(percentile = 0.05) |>
    posterior_percentile(percentile = 0.50, col_name = "eb50") |>
    posterior_percentile(percentile = 0.95) |>
    detect_signals(threshold = threshold) |>
    arrange(desc(eb05))
}

#' Filter results to signals only
#' @param results Output from run_signal_detection().
#' @return Filtered tibble.
#' @export
get_signals <- function(results) {
  results |>
    filter(.data$is_signal) |>
    select("drug", "event", "observed", "expected", "rr",
           "eb05", "eb50", "eb95", "signal_strength")
}
