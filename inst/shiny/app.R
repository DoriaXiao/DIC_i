# ===========================================================================
# DIC_i Interactive Demonstration — Factor Analysis
# v5: ggplot2 version for shinyapps.io deployment
# ===========================================================================

library(shiny)
library(ggplot2)

# ---- Single-chain Gibbs sampler ----
run_one_chain <- function(X, N, p, n_iter, burnin,
                          init_sign = 1, tau2 = 10,
                          a0 = 2, b0 = 1) {
  lam <- rep(init_sign * 0.5, p)
  psi <- rep(1.0, p)
  f   <- rnorm(N)
  n_keep <- n_iter - burnin
  lam_store <- matrix(NA_real_, n_keep, p)
  psi_store <- matrix(NA_real_, n_keep, p)
  ll_mat    <- matrix(NA_real_, n_keep, N)

  for (it in seq_len(n_iter)) {
    prec_f <- 1 + sum(lam^2 / psi)
    var_f  <- 1 / prec_f
    mean_f <- as.numeric(var_f * X %*% (lam / psi))
    f <- rnorm(N, mean_f, sqrt(var_f))
    ssf <- sum(f^2)
    for (j in seq_len(p)) {
      v <- 1 / (ssf / psi[j] + 1 / tau2)
      m <- v * sum(f * X[, j]) / psi[j]
      lam[j] <- rnorm(1, m, sqrt(v))
    }
    for (j in seq_len(p)) {
      r <- X[, j] - lam[j] * f
      psi[j] <- 1 / rgamma(1, a0 + N / 2, b0 + 0.5 * sum(r^2))
    }
    if (it > burnin) {
      k <- it - burnin
      lam_store[k, ] <- lam
      psi_store[k, ] <- psi
      Sigma <- tcrossprod(lam) + diag(psi, nrow = p)
      R <- tryCatch(chol(Sigma), error = function(e) NULL)
      if (is.null(R)) { ll_mat[k, ] <- NA; next }
      logdet <- 2 * sum(log(diag(R)))
      Si <- chol2inv(R)
      for (i in seq_len(N)) {
        ll_mat[k, i] <- -0.5 * (p * log(2 * pi) + logdet +
                                  sum(X[i, ] * (Si %*% X[i, ])))
      }
    }
  }
  list(lam = lam_store, psi = psi_store, ll = ll_mat)
}

# ---- Run 4 chains, combine, compute everything ----
run_multi_chain <- function(N, p, lambda_val, psi_val = 1.0,
                            n_iter = 600, burnin = 200, n_chains = 4) {
  lam_true <- rep(lambda_val, p)
  Sigma_true <- tcrossprod(lam_true) + diag(psi_val, p)
  R_true <- chol(Sigma_true)
  X <- matrix(rnorm(N * p), N, p) %*% R_true
  X <- scale(X, center = TRUE, scale = FALSE)

  signs <- c(1, 1, -1, -1)
  chains <- vector("list", n_chains)
  for (ch in seq_len(n_chains))
    chains[[ch]] <- run_one_chain(X, N, p, n_iter, burnin, init_sign = signs[ch])

  lam_all <- do.call(rbind, lapply(chains, "[[", "lam"))
  psi_all <- do.call(rbind, lapply(chains, "[[", "psi"))
  ll_all  <- do.call(rbind, lapply(chains, "[[", "ll"))
  ok <- complete.cases(ll_all)
  if (sum(ok) < 50) return(NULL)
  lam_all <- lam_all[ok, , drop = FALSE]
  psi_all <- psi_all[ok, , drop = FALSE]
  ll_all  <- ll_all[ok, , drop = FALSE]

  D_draws <- -2 * rowSums(ll_all)
  lam_bar <- colMeans(lam_all); psi_bar <- colMeans(psi_all)
  Sig_bar <- tcrossprod(lam_bar) + diag(psi_bar, nrow = p)
  R_bar   <- tryCatch(chol(Sig_bar), error = function(e) NULL)
  if (is.null(R_bar)) return(NULL)
  ld_bar  <- 2 * sum(log(diag(R_bar))); Si_bar <- chol2inv(R_bar)
  D_theta_bar <- -2 * sum(vapply(seq_len(N), function(i)
    -0.5 * (p*log(2*pi) + ld_bar + sum(X[i,] * (Si_bar %*% X[i,]))), numeric(1)))

  E_D       <- mean(D_draws)
  pD_val    <- E_D - D_theta_bar
  pV_val    <- 0.5 * var(D_draws)
  pwaic_val <- sum(apply(ll_all, 2, var))
  lppd <- sum(vapply(seq_len(N), function(i) {
    ll <- ll_all[,i]; mx <- max(ll); mx + log(mean(exp(ll - mx)))
  }, numeric(1)))

  data.frame(
    pD = pD_val, pV = pV_val, p_waic = pwaic_val,
    DIC    = D_theta_bar + 2 * pD_val,
    DIC_i  = E_D + pV_val,
    DIC_p  = D_theta_bar + 2 * pV_val,
    WAIC   = -2 * lppd + 2 * pwaic_val,
    k = 2L * p, N = N, p = p, lambda = lambda_val,
    sign_switched = {
      cs <- sapply(chains, function(ch) sign(mean(ch$lam[,1], na.rm=TRUE)))
      (sum(cs>0) > 0 & sum(cs<0) > 0)
    }
  )
}


# ==== Plotting helpers ====

# Plot 1: Penalty breakdown — faceted with free_y
plot_penalties <- function(df, x_var, x_lab) {
  k_true <- df$k[1]
  pen <- data.frame(
    xval = rep(df[[x_var]], 3),
    value = c(df$pV, df$p_waic, df$pD),
    Penalty = factor(
      rep(c("p_V (variance-based)", "p_WAIC", "p_DIC (plug-in)"), each = nrow(df)),
      levels = c("p_V (variance-based)", "p_WAIC", "p_DIC (plug-in)"))
  )
  ggplot(pen, aes(x = xval, y = value, colour = Penalty)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(se = FALSE, method = "loess", span = 0.55, linewidth = 1) +
    geom_hline(yintercept = k_true, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    facet_wrap(~ Penalty, ncol = 3, scales = "free_y") +
    scale_colour_manual(values = c(
      "p_V (variance-based)" = "#377EB8",
      "p_WAIC" = "#4DAF4A",
      "p_DIC (plug-in)" = "#E41A1C")) +
    labs(x = x_lab, y = "Effective Number of Parameters",
         title = "Penalty Stability: p_V and p_WAIC near k; p_DIC collapses") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 11),
          strip.background = element_rect(fill = "grey97"),
          plot.title = element_text(face = "bold", size = 14))
}

# Plot 2a: p_V and p_WAIC overlaid
plot_penalty_overlay <- function(df, x_var, x_lab) {
  k_true <- df$k[1]
  dat <- data.frame(
    xval = rep(df[[x_var]], 2),
    value = c(df$pV, df$p_waic),
    Penalty = factor(rep(c("p_V", "p_WAIC"), each = nrow(df)))
  )
  ggplot(dat, aes(x = xval, y = value, colour = Penalty, shape = Penalty)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(se = FALSE, method = "loess", span = 0.55, linewidth = 1) +
    geom_hline(yintercept = k_true, linetype = "dashed", colour = "black") +
    annotate("text", x = min(dat$xval), y = k_true + 0.3,
             label = paste0("k = ", k_true), hjust = 0, fontface = "italic", size = 4) +
    scale_colour_manual(values = c("p_V" = "#377EB8", "p_WAIC" = "#4DAF4A")) +
    scale_shape_manual(values = c(16, 17)) +
    labs(x = x_lab, y = "Effective Number of Parameters",
         title = "Both penalties converge to k",
         colour = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = c(0.85, 0.15),
          legend.background = element_rect(fill = alpha("white", 0.8)),
          plot.title = element_text(face = "bold", size = 13))
}

# Plot 2b: p_V - p_WAIC difference
plot_penalty_gap <- function(df, x_var, x_lab) {
  dat <- data.frame(xval = df[[x_var]], diff = df$pV - df$p_waic)
  ggplot(dat, aes(x = xval, y = diff)) +
    geom_point(colour = "#984EA3", alpha = 0.7, size = 2.5) +
    geom_smooth(se = FALSE, method = "loess", span = 0.55,
                colour = "#984EA3", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    labs(x = x_lab, y = expression(p[V] - p[WAIC]),
         title = expression(bold("Penalty gap:") ~ p[V] - p[WAIC] %->% 0)) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 13))
}

# Plot 3: Criterion alignment — delta from WAIC
plot_criterion_deltas <- function(df, x_var, x_lab) {
  dat <- data.frame(
    xval = rep(df[[x_var]], 3),
    delta = c(df$DIC - df$WAIC,
              df$DIC_p - df$WAIC,
              df$DIC_i - df$WAIC),
    Criterion = factor(
      rep(c("DIC (classical)", "DIC_p (Gelman)", "DIC_i (proposed)"),
          each = nrow(df)),
      levels = c("DIC (classical)", "DIC_p (Gelman)", "DIC_i (proposed)"))
  )
  ggplot(dat, aes(x = xval, y = delta, colour = Criterion)) +
    geom_point(alpha = 0.65, size = 2.5) +
    geom_smooth(se = FALSE, method = "loess", span = 0.55, linewidth = 1.1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.6) +
    scale_colour_manual(
      values = c("DIC (classical)" = "#E41A1C",
                 "DIC_p (Gelman)" = "#FF7F00",
                 "DIC_i (proposed)" = "#377EB8")) +
    labs(x = x_lab,
         y = "Criterion - WAIC",
         title = "Alignment with WAIC: DIC_i tracks WAIC; DIC and DIC_p diverge as mirror images",
         colour = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", size = 12))
}

# Summary table
make_summary <- function(df) {
  v <- list(
    p_DIC = df$pD, p_V = df$pV, p_WAIC = df$p_waic,
    DIC = df$DIC, DIC_p = df$DIC_p,
    DIC_i = df$DIC_i, WAIC = df$WAIC)
  data.frame(
    Quantity = names(v),
    Mean = sapply(v, mean), SD = sapply(v, sd),
    Min = sapply(v, min), Max = sapply(v, max),
    row.names = NULL)
}


# ====================================================================
# UI
# ====================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Palatino Linotype', 'Book Antiqua', Palatino, Georgia, serif;
           color: #2c3e50; }
    h2, h3, h4 { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }
    .btn-primary { background-color: #2c3e50; border-color: #1a252f; }
    .btn-primary:hover { background-color: #34495e; }
    .help-text { color: #7f8c8d; font-size: 0.85em; margin-top: 6px; line-height: 1.5; }
    .result-note { font-size: 0.92em; color: #2c3e50; margin-top: 10px;
                   padding: 10px 14px; background: #ecf0f1; border-left: 4px solid #3498db;
                   border-radius: 2px; }
  "))),

  titlePanel(HTML("DIC<sub><i>i</i></sub>: A Parameterization-Invariant DIC for Latent Variable Models")),
  p(HTML("Interactive companion to Xiao &amp; Rabe-Hesketh (2026), <a href='https://arxiv.org/abs/2605.27844' target='_blank'>arXiv:2605.27844</a>. Factor analysis under reflection invariance."),
    style = "color: #7f8c8d; margin-bottom: 12px;"),

  tabsetPanel(id = "main_tabs", type = "pills",

              # ===================== TAB 1 =====================
              tabPanel(HTML("Sensitivity to &lambda;"),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      sliderInput("n_tab1", "Sample size (N):", min = 200, max = 1000, value = 400, step = 50),
                                      sliderInput("p_tab1", "Items (p):", min = 4, max = 12, value = 6, step = 1),
                                      sliderInput("npts_tab1", HTML("Grid points:"), min = 10, max = 40, value = 20, step = 5),
                                      actionButton("run_tab1", "Generate New Data", class = "btn-primary"),
                                      hr(),
                                      p(class = "help-text", HTML(
                                        "<b>Design:</b> 1-factor, symmetric priors, 4 chains at opposite signs. The model has <i>k</i> = 2<i>p</i> free parameters (<i>p</i> loadings, <i>p</i> unique variances; data are mean-centered, so no intercepts).<br><br>
             <b>Mechanism:</b> When chains settle in different sign modes, E[&lambda;] is pulled toward zero &rarr;
             D(&theta;&#772;) explodes &rarr; p<sub>DIC</sub> goes massively negative.<br><br>
             <b>Key:</b> p<sub><i>V</i></sub> depends on &lambda;&lambda;&prime; (sign-invariant), stays near <i>k</i>."))
                         ),
                         mainPanel(width = 9,
                                   plotOutput("p1_penalties", height = "280px"),
                                   fluidRow(
                                     column(6, plotOutput("p1_overlay", height = "270px")),
                                     column(6, plotOutput("p1_gap", height = "270px"))
                                   ),
                                   plotOutput("p1_deltas", height = "300px"),
                                   hr(),
                                   h4("Summary"), tableOutput("t1_table"),
                                   div(class = "result-note", textOutput("t1_note"))
                         )
                       )
              ),

              # ===================== TAB 2 =====================
              tabPanel("Sensitivity to N",
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      sliderInput("lam_tab2", HTML("Loading (&lambda;):"), min = 0.3, max = 1.2, value = 0.7, step = 0.05),
                                      sliderInput("p_tab2", "Items (p):", min = 4, max = 12, value = 6, step = 1),
                                      sliderInput("npts_tab2", "Grid points:", min = 10, max = 30, value = 15, step = 5),
                                      actionButton("run_tab2", "Generate New Data", class = "btn-primary"),
                                      hr(),
                                      p(class = "help-text", HTML(
                                        "<b>Expect:</b> Larger N sharpens modes &rarr; p<sub>DIC</sub> gets worse.<br>
             p<sub><i>V</i></sub>, p<sub>WAIC</sub> converge to <i>k</i>.<br>
             DIC<sub><i>i</i></sub> &minus; WAIC &rarr; 0."))
                         ),
                         mainPanel(width = 9,
                                   plotOutput("p2_penalties", height = "280px"),
                                   fluidRow(
                                     column(6, plotOutput("p2_overlay", height = "270px")),
                                     column(6, plotOutput("p2_gap", height = "270px"))
                                   ),
                                   plotOutput("p2_deltas", height = "300px"),
                                   hr(),
                                   h4("Summary"), tableOutput("t2_table"),
                                   div(class = "result-note", textOutput("t2_note"))
                         )
                       )
              )
  )
)


# ====================================================================
# Server
# ====================================================================
server <- function(input, output, session) {

  res1 <- reactiveVal(NULL)
  res2 <- reactiveVal(NULL)

  # ---- Tab 1 ----
  observeEvent(input$run_tab1, {
    N_fix <- input$n_tab1; p_fix <- input$p_tab1; n_pts <- input$npts_tab1
    withProgress(message = "Running 4-chain Gibbs samplers...", value = 0, {
      lambdas <- sort(runif(n_pts, 0.2, 1.2))
      out <- vector("list", n_pts)
      for (i in seq_along(lambdas)) {
        incProgress(1/n_pts, detail = paste0(round(lambdas[i],2)," (",i,"/",n_pts,")"))
        out[[i]] <- run_multi_chain(N=N_fix, p=p_fix, lambda_val=lambdas[i])
      }
    })
    res1(do.call(rbind, Filter(Negate(is.null), out)))
  })

  output$p1_penalties <- renderPlot({ df <- res1(); if (is.null(df)) return(NULL)
  plot_penalties(df, "lambda", expression("Loading scale " * (lambda))) })
  output$p1_overlay   <- renderPlot({ df <- res1(); if (is.null(df)) return(NULL)
  plot_penalty_overlay(df, "lambda", expression("Loading scale " * (lambda))) })
  output$p1_gap       <- renderPlot({ df <- res1(); if (is.null(df)) return(NULL)
  plot_penalty_gap(df, "lambda", expression("Loading scale " * (lambda))) })
  output$p1_deltas    <- renderPlot({ df <- res1(); if (is.null(df)) return(NULL)
  plot_criterion_deltas(df, "lambda", expression("Loading scale " * (lambda))) })
  output$t1_table     <- renderTable({ df <- res1(); if (is.null(df)) return(NULL)
  make_summary(df) }, digits = 2)
  output$t1_note      <- renderText({
    df <- res1(); if (is.null(df)) return("Press 'Generate New Data' to begin.")
    n_sw <- sum(df$sign_switched)
    paste0("Sign switching: ", n_sw, "/", nrow(df), " (", round(100*n_sw/nrow(df)), "%). ",
           "Mean p_V = ", round(mean(df$pV),1), ", p_WAIC = ", round(mean(df$p_waic),1),
           " (k = ", df$k[1], "). |DIC_i - WAIC| = ", round(mean(abs(df$DIC_i-df$WAIC)),1), ".")
  })

  # ---- Tab 2 ----
  observeEvent(input$run_tab2, {
    lam_fix <- input$lam_tab2; p_fix <- input$p_tab2; n_pts <- input$npts_tab2
    withProgress(message = "Running 4-chain Gibbs samplers...", value = 0, {
      Ns <- sort(round(runif(n_pts, 200, 1000)))
      out <- vector("list", n_pts)
      for (i in seq_along(Ns)) {
        incProgress(1/n_pts, detail = paste0("N=", Ns[i], " (", i, "/", n_pts, ")"))
        out[[i]] <- run_multi_chain(N=Ns[i], p=p_fix, lambda_val=lam_fix)
      }
    })
    res2(do.call(rbind, Filter(Negate(is.null), out)))
  })

  output$p2_penalties <- renderPlot({ df <- res2(); if (is.null(df)) return(NULL)
  plot_penalties(df, "N", "Sample size (N)") })
  output$p2_overlay   <- renderPlot({ df <- res2(); if (is.null(df)) return(NULL)
  plot_penalty_overlay(df, "N", "Sample size (N)") })
  output$p2_gap       <- renderPlot({ df <- res2(); if (is.null(df)) return(NULL)
  plot_penalty_gap(df, "N", "Sample size (N)") })
  output$p2_deltas    <- renderPlot({ df <- res2(); if (is.null(df)) return(NULL)
  plot_criterion_deltas(df, "N", "Sample size (N)") })
  output$t2_table     <- renderTable({ df <- res2(); if (is.null(df)) return(NULL)
  make_summary(df) }, digits = 2)
  output$t2_note      <- renderText({
    df <- res2(); if (is.null(df)) return("Press 'Generate New Data' to begin.")
    n_sw <- sum(df$sign_switched)
    paste0("Sign switching: ", n_sw, "/", nrow(df), " (", round(100*n_sw/nrow(df)), "%). ",
           "Mean p_V = ", round(mean(df$pV),1), ", p_WAIC = ", round(mean(df$p_waic),1),
           " (k = ", df$k[1], "). |DIC_i - WAIC| = ", round(mean(abs(df$DIC_i-df$WAIC)),1), ".")
  })
}

shinyApp(ui = ui, server = server)
