#' DICi: Parameterization-Invariant Deviance Information Criterion
#'
#' Computes DIC_i, a parameterization-invariant, plug-in-free version of the
#' Deviance Information Criterion for Bayesian model comparison in latent
#' variable models (Xiao and Rabe-Hesketh, 2026).
#'
#' @section Main functions:
#' \describe{
#'   \item{[compute_dic_i()]}{Compute DIC_i from a log-likelihood matrix
#'     or deviance draws.}
#'   \item{[dic_i_from_cmdstanr()]}{Convenience wrapper for CmdStanR fit
#'     objects.}
#'   \item{[compare_dic_i()]}{Compare DIC_i across multiple models.}
#'   \item{[dici_example()]}{Access bundled example files (Stan model and
#'     R script).}
#' }
#'
#' @section Vignettes:
#' \describe{
#'   \item{`getting-started`}{Quick start with no dependencies beyond base R.}
#'   \item{`stan-workflow`}{Full factor analysis example in Stan with sign
#'     switching demonstration.}
#' }
#'
#' @references
#' Xiao, X. and Rabe-Hesketh, S. (2026). A Parameterization-Invariant DIC.
#'   *arXiv preprint* arXiv:2605.27844.
#'
#' @docType package
#' @name DICi-package
"_PACKAGE"
