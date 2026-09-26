# Extract DIC_i from a CmdStanR fit object

Convenience function that extracts the pointwise marginal log-likelihood
matrix from a cmdstanr fit object and computes \\\mathrm{DIC}\_i\\.

## Usage

``` r
dic_i_from_cmdstanr(fit, log_lik_name = "log_lik")
```

## Arguments

- fit:

  A `CmdStanMCMC` object, as returned by the `$sample()` method of a
  cmdstanr `CmdStanModel`.

- log_lik_name:

  Character string giving the name of the log-likelihood variable in the
  Stan model's `generated quantities` block. Default `"log_lik"`.

## Value

A `DICi` object (see
[`compute_dic_i()`](https://doriaxiao.github.io/DIC_i/reference/compute_dic_i.md)).

## Details

This function requires that:

1.  The `cmdstanr` and `posterior` packages are installed.

2.  The Stan model includes a `generated quantities` block that computes
    **marginal** pointwise log-likelihoods stored in an array named
    `log_lik_name`.

The function extracts the log-likelihood draws as an S x N matrix using
`fit$draws()` and passes it to
[`compute_dic_i()`](https://doriaxiao.github.io/DIC_i/reference/compute_dic_i.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(cmdstanr)
fit <- mod$sample(data = stan_data, chains = 4, iter_sampling = 1000)
dic_i_from_cmdstanr(fit)
} # }
```
