#' Extract DIC_i from a CmdStanR fit object
#'
#' Convenience function that extracts the pointwise marginal log-likelihood
#' matrix from a [cmdstanr] fit object and computes DIC_i.
#'
#' @param fit A `CmdStanMCMC` object from [cmdstanr::sample()].
#' @param log_lik_name Character string giving the name of the log-likelihood
#'   variable in the Stan model's `generated quantities` block. Default
#'   `"log_lik"`.
#'
#' @return A `DICi` object (see [compute_dic_i()]).
#'
#' @details
#' This function requires that:
#' 1. The `cmdstanr` and `posterior` packages are installed.
#' 2. The Stan model includes a `generated quantities` block that computes
#'    **marginal** pointwise log-likelihoods stored in an array named
#'    `log_lik_name`.
#'
#' The function extracts the log-likelihood draws as an S x N matrix using
#' `fit$draws()` and passes it to [compute_dic_i()].
#'
#' @examples
#' \dontrun{
#' library(cmdstanr)
#' fit <- mod$sample(data = stan_data, chains = 4, iter_sampling = 1000)
#' dic_i_from_cmdstanr(fit)
#' }
#'
#' @export
dic_i_from_cmdstanr <- function(fit, log_lik_name = "log_lik") {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required. Install from https://mc-stan.org/cmdstanr/")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required. Install with install.packages('posterior').")
  }

  log_lik <- fit$draws(log_lik_name, format = "draws_matrix")
  log_lik <- unclass(log_lik)
  attr(log_lik, "dimnames") <- NULL

  compute_dic_i(log_lik = log_lik)
}
