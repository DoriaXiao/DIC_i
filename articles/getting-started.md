# Getting Started with DICi

## Installation

``` r

remotes::install_github("DoriaXiao/DIC_i")
```

## What is DIC_i?

DIC_i is a parameterization-invariant, plug-in-free version of the
Deviance Information Criterion for Bayesian model comparison. It adds
the variance-based penalty of Gelman et al. (2014) to the posterior
**mean** deviance, rather than to the unstable plug-in deviance used by
the classic DIC:

\\\mathrm{DIC}\_i = E\[D(\theta)\] + \frac{1}{2}
\mathrm{Var}\[D(\theta)\]\\

where \\D(\theta) = -2 \log f_m(y \mid \theta)\\ is the marginal
deviance. The penalty \\p_V = \frac{1}{2}\mathrm{Var}(D)\\ is always
non-negative, and because DIC_i never evaluates the deviance at a point
estimate, the whole criterion is invariant to reparameterization.

This is what distinguishes DIC_i from the related variance-based DIC of
Gelman et al. (2014), \\\mathrm{DIC}\_p = D(\bar\theta) + 2 p_V\\, which
keeps the plug-in deviance \\D(\bar\theta)\\ and so remains unstable
under multimodality.

## Quick start

The core function takes an S x N matrix of pointwise marginal
log-likelihoods (rows = posterior draws, columns = observations):

``` r

library(DICi)

# Simulate a log-likelihood matrix (S=1000 draws, N=50 observations)
set.seed(42)
log_lik <- matrix(rnorm(1000 * 50, mean = -2, sd = 0.5), nrow = 1000, ncol = 50)

result <- compute_dic_i(log_lik)
print(result)
#> Parameterization-Invariant DIC (DIC_i)
#> --------------------------------------
#>   DIC_i  = 225.69
#>   p_V    = 25.55
#>   E[D]   = 200.14
#> --------------------------------------
```

## Using deviance draws directly

If you only have the joint marginal deviance (e.g., from Mplus output),
pass it directly:

``` r

deviance_draws <- -2 * rowSums(log_lik)
result2 <- compute_dic_i(deviance_draws = deviance_draws)
print(result2)
#> Parameterization-Invariant DIC (DIC_i)
#> --------------------------------------
#>   DIC_i  = 225.69
#>   p_V    = 25.55
#>   E[D]   = 200.14
#> --------------------------------------
```

The results are identical because DIC_i only needs \\E\[D\]\\ and
\\\mathrm{Var}(D)\\.

## Comparing models

Fit multiple candidate models, compute DIC_i for each, and compare:

``` r

set.seed(123)
ll_simple  <- matrix(rnorm(1000 * 50, mean = -2.0, sd = 0.4), 1000, 50)
ll_complex <- matrix(rnorm(1000 * 50, mean = -1.9, sd = 0.6), 1000, 50)

fit_simple  <- compute_dic_i(ll_simple)
fit_complex <- compute_dic_i(ll_complex)

compare_dic_i(simple = fit_simple, complex = fit_complex)
#>     model    dic_i      p_v      e_d delta_dic_i
#> 1  simple 215.8741 15.80440 200.0697     0.00000
#> 2 complex 227.7996 38.02143 189.7782    11.92546
```

The model with the lowest DIC_i is preferred. The `delta_dic_i` column
shows the difference from the best model.

## Important: Marginal log-likelihoods required

For latent variable models (factor analysis, mixture models, multilevel
models), the input must be **marginal** log-likelihoods, integrated over
the latent variables. Conditional log-likelihoods (evaluated at specific
latent variable draws) will produce misleading results.

See
[`vignette("stan-workflow")`](https://doriaxiao.github.io/DIC_i/articles/stan-workflow.md)
for a complete example of fitting a factor analysis model in Stan and
extracting marginal log-likelihoods.

## Next steps

- [`vignette("stan-workflow")`](https://doriaxiao.github.io/DIC_i/articles/stan-workflow.md)
  — Full Stan example with sign switching demonstration
- [`?compute_dic_i`](https://doriaxiao.github.io/DIC_i/reference/compute_dic_i.md)
  — Function reference with details on marginal vs. conditional
  likelihoods
- [`dici_example()`](https://doriaxiao.github.io/DIC_i/reference/dici_example.md)
  — Copy example Stan model and R script to your working directory
