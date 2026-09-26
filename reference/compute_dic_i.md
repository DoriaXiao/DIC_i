# Compute the Parameterization-Invariant DIC (DIC_i)

Computes DIC_i, the parameterization-invariant Deviance Information
Criterion proposed by Xiao and Rabe-Hesketh (2026). Unlike the classic
DIC, DIC_i does not depend on a plug-in estimate of the deviance and is
therefore robust to the multimodal posteriors caused by sign switching,
label switching, and parameterization switching in latent variable
models.

## Usage

``` r
compute_dic_i(log_lik = NULL, deviance_draws = NULL)
```

## Arguments

- log_lik:

  An S x N matrix of pointwise **marginal** log-likelihoods, where S is
  the number of posterior MCMC draws (rows) and N is the number of
  observations or clusters (columns). Each entry `log_lik[s, i]`
  represents \\\log f(y_i \mid \theta^{(s)})\\, the marginal
  log-likelihood contribution of observation \\i\\ evaluated at
  posterior draw \\s\\.

  **Important:** For latent variable models (e.g., factor analysis,
  mixed models, mixture models), the log-likelihoods must be marginal
  over the latent variables, not conditional on them. See Details.

- deviance_draws:

  Optional numeric vector of length S containing pre-computed marginal
  deviance draws \\D(\theta^{(s)}) = -2 \sum_i \log f(y_i \mid
  \theta^{(s)})\\. If provided, `log_lik` is ignored. This is useful
  when only the joint deviance is available (e.g., from Mplus output)
  and pointwise log-likelihoods are not.

## Value

A named list of class `"DICi"` with components:

- `dic_i`:

  The DIC_i value: \\\bar{D} + p_V\\.

- `p_v`:

  The variance-based effective number of parameters \\p_V = \frac{1}{2}
  \mathrm{Var}(D(\theta))\\ of Gelman et al. (2014).

- `e_d`:

  The posterior mean deviance \\\bar{D} = E\[D(\theta)\]\\.

## Details

DIC_i is defined as: \$\$\mathrm{DIC}\_i = E\_{\theta\|y}\[D(\theta)\] +
\frac{1}{2} \mathrm{Var}\_{\theta\|y}\[D(\theta)\]\$\$

where \\D(\theta) = -2 \log f_m(y \mid \theta)\\ is the marginal
deviance. The criterion adds the variance-based penalty \\p_V =
\frac{1}{2} \mathrm{Var}(D)\\ of Gelman et al. (2014) to the posterior
mean deviance. Both terms are invariant to reparameterization of the
model parameters, and \\p_V\\ is always non-negative.

The "i" in DIC_i stands for (parameterization-)**i**nvariant. The key
difference from the related variance-based DIC of Gelman et al. (2014),
here denoted \\\mathrm{DIC}\_p = D(\bar{\theta}) + 2 p_V\\, is that
DIC_i uses the posterior **mean** deviance \\E\[D(\theta)\]\\ in place
of the plug-in deviance \\D(\bar{\theta})\\. The plug-in deviance is
what makes the classic DIC (and DIC_p) unstable under multimodality, so
removing it is what gives DIC_i its invariance. DIC_i is asymptotically
equivalent to the WAIC (Watanabe, 2010) but does not require the
likelihood to factorize into independent pointwise contributions.

### Marginal vs. conditional log-likelihoods

For meaningful results in latent variable models, the input must be
**marginal** log-likelihoods (integrated over latent variables), not
conditional log-likelihoods (evaluated at specific latent variable
draws).

- **Stan:** If you code the marginal likelihood in the `model` block
  (e.g., using `log_sum_exp` for mixtures or the multivariate normal
  density after integrating out random effects) and output pointwise
  contributions in `generated quantities`, the extracted log-likelihoods
  are marginal. This is the recommended workflow.

- **JAGS/BUGS:** These programs typically sample latent variables and
  report conditional likelihoods. To use DIC_i, you would need to
  compute the marginal likelihood yourself, which requires
  model-specific integration (analytic for linear mixed models,
  numerical otherwise).

- **Mplus:** The marginal deviance is reported directly. Use the
  `deviance_draws` argument to pass in the deviance values.

For further discussion of marginal vs. conditional likelihoods in
Bayesian model comparison, see Merkle, Furr, and Rabe-Hesketh (2019).

## References

Xiao, X. and Rabe-Hesketh, S. (2026). A Parameterization-Invariant DIC.
*arXiv preprint* arXiv:2605.27844.

Gelman, A., Hwang, J., and Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models. *Statistics and Computing*,
24, 997–1016.

Merkle, E.C., Furr, D., and Rabe-Hesketh, S. (2019). Bayesian comparison
of latent variable models: Conditional versus marginal likelihoods.
*Psychometrika*, 84, 802–829.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research*, 11, 3571–3594.

## Examples

``` r
# Simulate fake log-likelihood matrix (S=1000 draws, N=50 observations)
set.seed(42)
S <- 1000
N <- 50
log_lik <- matrix(rnorm(S * N, mean = -2, sd = 0.5), nrow = S, ncol = N)

result <- compute_dic_i(log_lik)
result$dic_i
#> [1] 225.6861
result$p_v
#> [1] 25.54797
result$e_d
#> [1] 200.1382

# Using pre-computed deviance draws
dev_draws <- -2 * rowSums(log_lik)
result2 <- compute_dic_i(deviance_draws = dev_draws)
all.equal(result$dic_i, result2$dic_i)
#> [1] TRUE
```
