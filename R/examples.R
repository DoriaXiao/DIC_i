#' Access Example Files
#'
#' Returns the path to example files bundled with the package, or copies
#' them to the current working directory.
#'
#' @param name Name of the example file. Available files:
#'   \describe{
#'     \item{`"fa_marginal.stan"`}{Stan model for one-factor analysis with
#'       marginal likelihood and pointwise log-likelihoods.}
#'     \item{`"example_fa.R"`}{Complete R script: simulate data, fit in
#'       Stan, compute DIC_i, compare with WAIC/LOO.}
#'   }
#'   If `NULL` (the default), lists all available example files.
#'
#' @param copy Logical. If `TRUE` (default), copies the file to the current
#'   working directory and returns the destination path. If `FALSE`, returns
#'   the path to the file inside the package without copying.
#'
#' @return The file path (invisibly when `copy = TRUE`).
#'
#' @examples
#' # List available examples
#' dici_example()
#'
#' # Get path without copying
#' dici_example("fa_marginal.stan", copy = FALSE)
#'
#' \dontrun{
#' # Copy to working directory (interactive use)
#' dici_example("fa_marginal.stan")
#' dici_example("example_fa.R")
#' }
#'
#' @export
dici_example <- function(name = NULL, copy = TRUE) {

  ex_dir <- system.file("examples", package = "DICi")

  if (ex_dir == "") {
    stop("Example files not found. Is the package installed?")
  }

  available <- list.files(ex_dir)

  if (is.null(name)) {
    cat("Available example files:\n")
    for (f in available) {
      cat("  ", f, "\n")
    }
    cat("\nUse dici_example(\"filename\") to copy to your working directory.\n")
    return(invisible(available))
  }

  src <- file.path(ex_dir, name)

  if (!file.exists(src)) {
    stop(
      sprintf("'%s' not found. Available files: %s",
              name, paste(available, collapse = ", "))
    )
  }

  if (!copy) {
    return(src)
  }

  dst <- file.path(getwd(), name)

  if (file.exists(dst)) {
    message(sprintf("'%s' already exists in working directory. Overwriting.", name))
  }

  file.copy(src, dst, overwrite = TRUE)
  message(sprintf("Copied '%s' to %s", name, dst))
  invisible(dst)
}
