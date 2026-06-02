# DIC_i interactive demo (Shiny)

Interactive companion to Xiao & Rabe-Hesketh (2026),
"A Parameterization-Invariant DIC" ([arXiv:2605.27844](https://arxiv.org/abs/2605.27844)).

Live: https://doriaxiao.shinyapps.io/dicv_app/

The app runs live 4-chain Gibbs samplers for a one-factor model (symmetric
priors, chains initialised at opposite signs to induce sign switching) and
shows, across a grid of loadings (Tab 1) or sample sizes (Tab 2):

1. **Penalty stability** — `p_V` and `p_WAIC` stay near the true `k`, while the
   classical plug-in penalty `p_DIC` collapses to large negative values.
2. **Penalty convergence** — `p_V - p_WAIC -> 0` (Section 4.2).
3. **Criterion alignment** — `DIC_i ≈ WAIC`, while classical `DIC` and Gelman's
   `DIC_p` diverge from WAIC as mirror images.

## Run locally

```r
# install.packages(c("shiny", "ggplot2"))
shiny::runApp(system.file("shiny", package = "DICi"))
# or, from a clone of this repo:
shiny::runApp("inst/shiny")
```

## Deploy (shinyapps.io)

```r
rsconnect::deployApp("inst/shiny", appName = "dicv_app", forceUpdate = TRUE)
```
