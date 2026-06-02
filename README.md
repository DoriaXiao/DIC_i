# DICi

Parameterization-Invariant Deviance Information Criterion (DIC_i) for Bayesian Model Comparison

## Overview

The `DICi` package computes **DIC_i**, a parameterization-invariant,
plug-in-free version of the Deviance Information Criterion for comparing
latent variable models (Xiao and Rabe-Hesketh, 2026,
[arXiv:2605.27844](https://arxiv.org/abs/2605.27844)).

DIC_i adds the variance-based penalty $p_V = \tfrac12 \mathrm{Var}(D)$ of
Gelman et al. (2014b) to the posterior **mean** deviance $E[D(\theta)]$:

$$\mathrm{DIC}_i = E[D(\theta)] + \tfrac12 \mathrm{Var}[D(\theta)].$$

Unlike the classical DIC, DIC_i never evaluates the deviance at a point
estimate, so it:

- **Never produces negative penalty terms**, even under sign switching
  (factor analysis), label switching (mixture models), or parameterization
  switching (overfitted mixtures)
- **Does not depend on plug-in estimates**, so it is not destabilized by
  multimodal posteriors
- **Is asymptotically equivalent to WAIC**, but does not require pointwise
  likelihood factorization
- **Is computationally simple**: only the posterior mean and variance of the
  marginal deviance are needed

### How DIC_i relates to other criteria

| Criterion | Formula | Source |
|---|---|---|
| **DIC_i** (this package) | $E[D(\theta)] + p_V$ | Xiao & Rabe-Hesketh (2026) |
| Classical DIC | $D(\bar\theta) + 2 p_{\mathrm{DIC}}$,  $p_{\mathrm{DIC}} = E[D] - D(\bar\theta)$ | Spiegelhalter et al. (2002) |
| DIC_p (variance-based) | $D(\bar\theta) + 2 p_V$ | Gelman et al. (2014b) |

> **Note.** The penalty $p_V = \tfrac12\mathrm{Var}(D)$ is **Gelman et al.
> (2014b)'s** penalty — it is not new here. The contribution of DIC_i is to
> drop the plug-in deviance $D(\bar\theta)$ and use the posterior mean
> deviance $E[D(\theta)]$ instead, which is what makes the criterion
> invariant to reparameterization. The closely related **DIC_p** of Gelman
> et al. (2014b) keeps $D(\bar\theta)$ and is therefore **not** robust to the
> multimodality this package targets — do not confuse the two.

## Installation

```r
remotes::install_github("DoriaXiao/DIC_i")
```

## Quick start

```r
library(DICi)

# From a log-likelihood matrix (S draws x N observations)
set.seed(42)
log_lik <- matrix(rnorm(1000 * 50, -2, 0.5), 1000, 50)
compute_dic_i(log_lik)
#> Parameterization-Invariant DIC (DIC_i)
#> --------------------------------------
#>   DIC_i  = 211.08
#>   p_V    = 10.43
#>   E[D]   = 200.65
#> --------------------------------------
```

## After fitting a Stan model

```r
library(cmdstanr)
fit <- mod$sample(data = stan_data, chains = 4, iter_sampling = 1000)

# One-line computation
result <- dic_i_from_cmdstanr(fit)
print(result)
```

## Comparing models

```r
fit_k1 <- mod_k1$sample(data = stan_data, ...)
fit_k2 <- mod_k2$sample(data = stan_data, ...)

compare_dic_i(
  K1 = dic_i_from_cmdstanr(fit_k1),
  K2 = dic_i_from_cmdstanr(fit_k2)
)
#>   model   dic_i   p_v     e_d delta_dic_i
#> 1    K2 2801.44 13.21 2788.23        0.00
#> 2    K1 2847.31 19.02 2828.29       45.87
```

## Side-by-side with WAIC and LOO-CV

```r
library(loo)
log_lik <- fit$draws("log_lik", format = "draws_matrix")

dic_result  <- compute_dic_i(unclass(log_lik))
waic_result <- waic(log_lik)
loo_result  <- loo(log_lik)

cat(sprintf("p_V = %.1f, p_WAIC = %.1f, p_LOO = %.1f\n",
            dic_result$p_v,
            waic_result$estimates["p_waic", "Estimate"],
            loo_result$estimates["p_loo", "Estimate"]))
# p_V = 19.0, p_WAIC = 17.4, p_LOO = 17.4
```

## Running the full example

The package includes a complete factor analysis example demonstrating sign
switching. The example simulates data, fits a one-factor model in Stan with
symmetric priors, computes DIC_i alongside classical DIC / DIC_p / WAIC / LOO,
and shows that DIC_i remains stable while the classical DIC penalty becomes
negative.

```r
library(DICi)

# List available examples
dici_example()

# Copy Stan model and R script to your working directory
dici_example("fa_marginal.stan")
dici_example("example_fa.R")

# Run the full demo (requires cmdstanr)
source("example_fa.R")
```

Expected output includes a comparison table:

```
==========================================
  Information Criteria Comparison
==========================================
  Plug-in-free criteria:
    DIC_i      =   5724.5  (p_V    =  17.7)
    WAIC       =   5725.0  (p_WAIC =  18.0)
    LOO-CV     =   5725.1  (p_LOO  =  18.0)
  ------------------------------------------
  Plug-in-dependent criteria:
    DIC (class)=   5706.8  (p_DIC  =  -0.0)
    DIC_p      =   5742.3   [Gelman et al., 2014b]
  ------------------------------------------
  True k       = 18
==========================================
```

## Marginal log-likelihoods required

For latent variable models, the input must be **marginal** log-likelihoods
(integrated over latent variables). Conditional log-likelihoods will produce
misleading results.

| Software | How to get marginal log-lik |
|---|---|
| **Stan** | Code the marginal likelihood in `model` block; output `log_lik` in `generated quantities` |
| **Mplus** | Marginal deviance is reported directly; pass to `compute_dic_i(deviance_draws = ...)` |
| **JAGS/BUGS** | Default output is conditional; marginal requires model-specific integration |

See `vignette("stan-workflow")` for details.

## Citation

```
Xiao, X. and Rabe-Hesketh, S. (2026). A Parameterization-Invariant DIC.
arXiv preprint arXiv:2605.27844.
```

## License

MIT
