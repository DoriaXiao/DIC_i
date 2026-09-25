# Access Example Files

Returns the path to example files bundled with the package, or copies
them to the current working directory.

## Usage

``` r
dici_example(name = NULL, copy = TRUE)
```

## Arguments

- name:

  Name of the example file. Available files:

  `"fa_marginal.stan"`

  :   Stan model for one-factor analysis with marginal likelihood and
      pointwise log-likelihoods.

  `"example_fa.R"`

  :   Complete R script: simulate data, fit in Stan, compute DIC_i,
      compare with WAIC/LOO.

  If `NULL` (the default), lists all available example files.

- copy:

  Logical. If `TRUE` (default), copies the file to the current working
  directory and returns the destination path. If `FALSE`, returns the
  path to the file inside the package without copying.

## Value

The file path (invisibly when `copy = TRUE`).

## Examples

``` r
# List available examples
dici_example()
#> Available example files:
#>    example_fa.R 
#>    fa_marginal.stan 
#> 
#> Use dici_example("filename") to copy to your working directory.

# Get path without copying
dici_example("fa_marginal.stan", copy = FALSE)
#> [1] "/home/runner/work/_temp/Library/DICi/examples/fa_marginal.stan"

if (FALSE) { # \dontrun{
# Copy to working directory (interactive use)
dici_example("fa_marginal.stan")
dici_example("example_fa.R")
} # }
```
