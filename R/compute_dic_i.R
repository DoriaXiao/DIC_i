#' Compute the Parameterization-Invariant DIC (DIC_i)
#'
#' Computes DIC_i, the parameterization-invariant Deviance Information
#' Criterion proposed by Xiao and Rabe-Hesketh (2026). Unlike the classical
#' DIC, DIC_i does not depend on a plug-in estimate of the deviance and is
#' therefore robust to the multimodal posteriors caused by sign switching,
#' label switching, and parameterization switching in latent variable models.
#'
#' @param log_lik An S x N matrix of pointwise **marginal** log-likelihoods,
#'   where S is the number of posterior MCMC draws (rows) and N is the number
#'   of observations or clusters (columns). Each entry `log_lik[s, i]`
#'   represents \eqn{\log f(y_i \mid \theta^{(s)})}, the marginal
#'   log-likelihood contribution of observation \eqn{i} evaluated at posterior
#'   draw \eqn{s}.
#'
#'   **Important:** For latent variable models (e.g., factor analysis, mixed
#'   models, mixture models), the log-likelihoods must be marginal over the
#'   latent variables, not conditional on them. See Details.
#'
#' @param deviance_draws Optional numeric vector of length S containing
#'   pre-computed marginal deviance draws \eqn{D(\theta^{(s)}) = -2 \sum_i
#'   \log f(y_i \mid \theta^{(s)})}. If provided, `log_lik` is ignored. This
#'   is useful when only the joint deviance is available (e.g., from Mplus
#'   output) and pointwise log-likelihoods are not.
#'
#' @return A named list of class `"DICi"` with components:
#'   \describe{
#'     \item{`dic_i`}{The DIC_i value: \eqn{\bar{D} + p_V}.}
#'     \item{`p_v`}{The variance-based effective number of parameters
#'       \eqn{p_V = \frac{1}{2} \mathrm{Var}(D(\theta))} of Gelman et al.
#'       (2014b).}
#'     \item{`e_d`}{The posterior mean deviance \eqn{\bar{D} =
#'       E[D(\theta)]}.}
#'   }
#'
#' @details
#' DIC_i is defined as:
#' \deqn{\mathrm{DIC}_i = E_{\theta|y}[D(\theta)] + \frac{1}{2}
#'   \mathrm{Var}_{\theta|y}[D(\theta)]}
#'
#' where \eqn{D(\theta) = -2 \log f_m(y \mid \theta)} is the marginal
#' deviance. The criterion adds the variance-based penalty
#' \eqn{p_V = \frac{1}{2} \mathrm{Var}(D)} of Gelman et al. (2014b) to the
#' posterior mean deviance. Both terms are always non-negative and invariant
#' to reparameterization of the model parameters.
#'
#' The "i" in DIC_i stands for (parameterization-)**i**nvariant. The key
#' difference from the related variance-based DIC of Gelman et al. (2014b),
#' here denoted \eqn{\mathrm{DIC}_p = D(\bar{\theta}) + 2 p_V}, is that DIC_i
#' uses the posterior **mean** deviance \eqn{E[D(\theta)]} in place of the
#' plug-in deviance \eqn{D(\bar{\theta})}. The plug-in deviance is what makes
#' the classical DIC (and DIC_p) unstable under multimodality, so removing it
#' is what gives DIC_i its invariance. DIC_i is asymptotically equivalent to
#' the WAIC (Watanabe, 2010) but does not require the likelihood to factorize
#' into independent pointwise contributions.
#'
#' ## Marginal vs. conditional log-likelihoods
#'
#' For meaningful results in latent variable models, the input must be
#' **marginal** log-likelihoods (integrated over latent variables), not
#' conditional log-likelihoods (evaluated at specific latent variable draws).
#'
#' - **Stan:** If you code the marginal likelihood in the `model` block
#'   (e.g., using `log_sum_exp` for mixtures or the multivariate normal
#'   density after integrating out random effects) and output pointwise
#'   contributions in `generated quantities`, the extracted log-likelihoods
#'   are marginal. This is the recommended workflow.
#'
#' - **JAGS/BUGS:** These programs typically sample latent variables and
#'   report conditional likelihoods. To use DIC_i, you would need to compute
#'   the marginal likelihood yourself, which requires model-specific
#'   integration (analytic for linear mixed models, numerical otherwise).
#'
#' - **Mplus:** The marginal deviance is reported directly. Use the
#'   `deviance_draws` argument to pass in the deviance values.
#'
#' For further discussion of marginal vs. conditional likelihoods in
#' Bayesian model comparison, see Merkle, Furr, and Rabe-Hesketh (2019).
#'
#' @references
#' Xiao, X. and Rabe-Hesketh, S. (2026). A Parameterization-Invariant DIC.
#'   *arXiv preprint* arXiv:2605.27844.
#'
#' Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive
#'   information criteria for Bayesian models. *Statistics and Computing*,
#'   24, 997--1016.
#'
#' Merkle, E.C., Furr, D., and Rabe-Hesketh, S. (2019). Bayesian comparison
#'   of latent variable models: Conditional versus marginal likelihoods.
#'   *Psychometrika*, 84, 802--829.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#'   widely applicable information criterion in singular learning theory.
#'   *Journal of Machine Learning Research*, 11, 3571--3594.
#'
#' @examples
#' # Simulate fake log-likelihood matrix (S=1000 draws, N=50 observations)
#' set.seed(42)
#' S <- 1000
#' N <- 50
#' log_lik <- matrix(rnorm(S * N, mean = -2, sd = 0.5), nrow = S, ncol = N)
#'
#' result <- compute_dic_i(log_lik)
#' result$dic_i
#' result$p_v
#' result$e_d
#'
#' # Using pre-computed deviance draws
#' dev_draws <- -2 * rowSums(log_lik)
#' result2 <- compute_dic_i(deviance_draws = dev_draws)
#' all.equal(result$dic_i, result2$dic_i)
#'
#' @importFrom stats var
#' @export
compute_dic_i <- function(log_lik = NULL, deviance_draws = NULL) {

  # --- Input validation ---
  if (is.null(log_lik) && is.null(deviance_draws)) {
    stop("Either 'log_lik' or 'deviance_draws' must be provided.")
  }

  if (!is.null(deviance_draws)) {
    if (!is.numeric(deviance_draws) || !is.vector(deviance_draws)) {
      stop("'deviance_draws' must be a numeric vector.")
    }
    if (length(deviance_draws) < 2) {
      stop("'deviance_draws' must have at least 2 elements.")
    }
    d <- deviance_draws
  } else {
    if (!is.matrix(log_lik) && !is.array(log_lik)) {
      stop("'log_lik' must be a matrix (S rows x N columns).")
    }
    log_lik <- as.matrix(log_lik)
    if (nrow(log_lik) < 2) {
      stop("'log_lik' must have at least 2 rows (posterior draws).")
    }
    d <- -2 * rowSums(log_lik)
  }

  # --- Compute DIC_i = E[D] + p_V, with p_V = 0.5 * Var(D) ---
  e_d   <- mean(d)
  p_v   <- 0.5 * var(d)
  dic_i <- e_d + p_v

  structure(
    list(
      dic_i = dic_i,
      p_v   = p_v,
      e_d   = e_d
    ),
    class = "DICi"
  )
}


#' Compare DIC_i Across Models
#'
#' Compares DIC_i values from multiple fitted models, returning a summary
#' table sorted by DIC_i (lowest = preferred).
#'
#' @param ... Named `DICi` objects (output of [compute_dic_i()]).
#'
#' @return A data frame with columns `model`, `dic_i`, `p_v`, `e_d`,
#'   `delta_dic_i` (difference from the best model).
#'
#' @examples
#' set.seed(42)
#' ll_1 <- matrix(rnorm(500 * 50, -2.0, 0.5), 500, 50)
#' ll_2 <- matrix(rnorm(500 * 50, -2.1, 0.5), 500, 50)
#'
#' fit1 <- compute_dic_i(ll_1)
#' fit2 <- compute_dic_i(ll_2)
#'
#' compare_dic_i(model_1 = fit1, model_2 = fit2)
#'
#' @export
compare_dic_i <- function(...) {
  models <- list(...)

  if (length(models) < 2) {
    stop("At least two models are required for comparison.")
  }

  nms <- names(models)
  if (is.null(nms) || any(nms == "")) {
    stop("All arguments must be named (e.g., compare_dic_i(K1 = fit1, K2 = fit2)).")
  }

  for (nm in nms) {
    if (!inherits(models[[nm]], "DICi")) {
      stop(sprintf("'%s' is not a DICi object. Use compute_dic_i() first.", nm))
    }
  }

  tbl <- data.frame(
    model    = nms,
    dic_i    = vapply(models, function(x) x$dic_i, numeric(1)),
    p_v      = vapply(models, function(x) x$p_v, numeric(1)),
    e_d      = vapply(models, function(x) x$e_d, numeric(1)),
    stringsAsFactors = FALSE
  )

  tbl <- tbl[order(tbl$dic_i), ]
  tbl$delta_dic_i <- tbl$dic_i - tbl$dic_i[1]
  rownames(tbl) <- NULL

  tbl
}


#' Print method for DICi objects
#'
#' @param x A `DICi` object.
#' @param digits Number of decimal places. Default 2.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.DICi <- function(x, digits = 2, ...) {
  cat("Parameterization-Invariant DIC (DIC_i)\n")
  cat("--------------------------------------\n")
  cat(sprintf("  DIC_i  = %.*f\n", digits, x$dic_i))
  cat(sprintf("  p_V    = %.*f\n", digits, x$p_v))
  cat(sprintf("  E[D]   = %.*f\n", digits, x$e_d))
  cat("--------------------------------------\n")
  invisible(x)
}
