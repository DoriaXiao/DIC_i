# ==========================================================================
#  DIC_i Example: One-Factor Model with Sign Switching
#  Full workflow: simulate -> fit in Stan -> compute DIC_i -> compare
# ==========================================================================

library(cmdstanr)
library(posterior)
library(loo)
library(DICi)

# --------------------------------------------------------------------------
# 1. Simulate data from a one-factor model
# --------------------------------------------------------------------------
set.seed(2026)
N <- 400      # sample size
P <- 6        # number of indicators

# True parameters
lambda_true <- c(0.81, 0.72, 0.63, 0.54, 0.45, 0.36)  # c=0.9 scaling
alpha_true  <- rep(0, P)
sigma_true  <- rep(sqrt(0.5), P)  # residual SD

# Generate data: y_j = alpha + lambda * eta_j + epsilon_j
eta <- rnorm(N)
Y <- matrix(NA, N, P)
for (j in 1:N) {
  Y[j, ] <- alpha_true + lambda_true * eta[j] + rnorm(P, 0, sigma_true)
}

cat("Data generated: N =", N, ", P =", P, "\n")
cat("True number of parameters: k =", 3 * P, "(6 intercepts + 6 loadings + 6 residual SDs)\n\n")

# --------------------------------------------------------------------------
# 2. Fit the model in Stan (marginal likelihood, symmetric priors)
# --------------------------------------------------------------------------
stan_data <- list(N = N, P = P, Y = Y)

mod <- cmdstan_model("fa_marginal.stan")

fit <- mod$sample(
  data            = stan_data,
  chains          = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  parallel_chains = 4,
  seed            = 42
)

# --------------------------------------------------------------------------
# 3. Check for sign switching
# --------------------------------------------------------------------------
cat("\n--- Posterior summary for loadings ---\n")
print(fit$summary("lambda"))

lambda_array <- fit$draws("lambda", format = "draws_array")

cat("\n--- R-hat diagnostics for lambda[1] ---\n")
lambda1    <- extract_variable_matrix(lambda_array, "lambda[1]")  # iterations x chains
rhat_rank  <- rhat(lambda1)
rhat_split <- rhat_basic(lambda1)  # basic split-R-hat (BDA3)
cat(sprintf("  Rank-normalized R-hat: %.3f\n", rhat_rank))
cat(sprintf("  Basic split R-hat:     %.3f\n", rhat_split))

if (rhat_split > 1.5) {
  cat("  >>> Large basic split R-hat indicates sign switching between chains.\n")
}

# --------------------------------------------------------------------------
# 4. Compute DIC_i
# --------------------------------------------------------------------------
result <- dic_i_from_cmdstanr(fit)
cat("\n")
print(result)

# --------------------------------------------------------------------------
# 5. Compute classic DIC and Gelman's DIC_p for comparison
# --------------------------------------------------------------------------
log_lik <- fit$draws("log_lik", format = "draws_matrix")
log_lik <- unclass(log_lik)

# Posterior mean deviance
D_bar <- mean(-2 * rowSums(log_lik))

# Plug-in deviance D(theta-bar): the marginal deviance evaluated at the
# posterior means of the parameters. (Averaging the pointwise log-likelihoods
# instead, -2 * sum(colMeans(log_lik)), equals D_bar exactly, so p_DIC would
# always be 0.)
alpha_bar  <- colMeans(fit$draws("alpha",  format = "draws_matrix"))
lambda_bar <- colMeans(fit$draws("lambda", format = "draws_matrix"))
sigma_bar  <- colMeans(fit$draws("sigma",  format = "draws_matrix"))
V_bar <- tcrossprod(lambda_bar) + diag(sigma_bar^2)
R_bar <- chol(V_bar)                                    # V_bar = t(R_bar) %*% R_bar
Z_bar <- backsolve(R_bar, t(Y) - alpha_bar, transpose = TRUE)
D_plugin <- N * P * log(2 * pi) + 2 * N * sum(log(diag(R_bar))) + sum(Z_bar^2)

# Classic penalty and DIC (Spiegelhalter et al., 2002): can go negative
p_DIC <- D_bar - D_plugin
DIC_classical <- D_plugin + 2 * p_DIC   # = D_bar + p_DIC

# DIC_p (Gelman et al., 2014): plug-in deviance + variance penalty.
# Distinct from DIC_i, which uses the posterior mean deviance instead.
DIC_p <- D_plugin + 2 * result$p_v

# --------------------------------------------------------------------------
# 6. Compute WAIC and LOO-CV
# --------------------------------------------------------------------------
waic_result <- waic(log_lik)
loo_result  <- loo(log_lik)

# --------------------------------------------------------------------------
# 7. Full comparison table
# --------------------------------------------------------------------------
cat("\n==========================================\n")
cat("  Information Criteria Comparison\n")
cat("==========================================\n")
cat("  Plug-in-free criteria:\n")
cat(sprintf("    DIC_i      = %8.1f  (p_V    = %5.1f)\n", result$dic_i, result$p_v))
cat(sprintf("    WAIC       = %8.1f  (p_WAIC = %5.1f)\n",
            waic_result$estimates["waic", "Estimate"],
            waic_result$estimates["p_waic", "Estimate"]))
cat(sprintf("    LOO-CV     = %8.1f  (p_LOO  = %5.1f)\n",
            loo_result$estimates["looic", "Estimate"],
            loo_result$estimates["p_loo", "Estimate"]))
cat("  ------------------------------------------\n")
cat("  Plug-in-dependent criteria:\n")
cat(sprintf("    DIC classic= %8.1f  (p_DIC  = %5.1f)\n", DIC_classical, p_DIC))
cat(sprintf("    DIC_p      = %8.1f   [Gelman et al., 2014]\n", DIC_p))
cat("  ------------------------------------------\n")
cat(sprintf("  True k       = %d\n", 3 * P))
cat("==========================================\n")

if (p_DIC < 0) {
  cat("\n  >>> p_DIC is NEGATIVE: sign switching detected!\n")
  cat("  >>> Classic DIC and Gelman's DIC_p are meaningless (both use the\n")
  cat("  >>> plug-in deviance D(theta-bar)).\n")
  cat("  >>> DIC_i remains stable and close to WAIC.\n")
} else {
  cat("\n  >>> No sign switching in this run.\n")
  cat("  >>> All criteria should agree closely.\n")
}

# --------------------------------------------------------------------------
# 8. Compare multiple models (template)
# --------------------------------------------------------------------------
# If you had a second model:
# fit_2f <- mod_2f$sample(data = stan_data, ...)
# dic_2f <- dic_i_from_cmdstanr(fit_2f)
# compare_dic_i(one_factor = result, two_factor = dic_2f)
