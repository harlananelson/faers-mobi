# faers.mobi — Vaccine Safety Signal Detection

## Overview
Rhino/Shiny app for Bayesian disproportionality analysis of vaccine adverse events using VAERS data. Uses the `safetysignal` R package as its statistical engine.

## Architecture
- **Framework:** Rhino (production Shiny)
- **Engine:** `safetysignal` package (2-component Gamma-Poisson)
- **Data source:** VAERS (Vaccine Adverse Event Reporting System)
- **Domain:** faers.mobi

## Key Files
- `app/main.R` — App entry point
- `app/logic/signal_engine.R` — Wraps safetysignal for this app
- `app/view/signal_table.R` — Signal results table module

## Related Projects
- `safetysignal` — Shared Bayesian engine package
- `aers-mobi` — Drug/device version
- `globalpatientsafety` — Broader safety dashboard
