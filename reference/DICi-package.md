# DICi: Parameterization-Invariant Deviance Information Criterion

Computes DIC_i, a parameterization-invariant, plug-in-free version of
the Deviance Information Criterion for Bayesian model comparison in
latent variable models (Xiao and Rabe-Hesketh, 2026).

## Main functions

- [`compute_dic_i()`](https://doriaxiao.github.io/DIC_i/reference/compute_dic_i.md):

  Compute DIC_i from a log-likelihood matrix or deviance draws.

- [`dic_i_from_cmdstanr()`](https://doriaxiao.github.io/DIC_i/reference/dic_i_from_cmdstanr.md):

  Convenience wrapper for CmdStanR fit objects.

- [`compare_dic_i()`](https://doriaxiao.github.io/DIC_i/reference/compare_dic_i.md):

  Compare DIC_i across multiple models.

- [`dici_example()`](https://doriaxiao.github.io/DIC_i/reference/dici_example.md):

  Access bundled example files (Stan model and R script).

## Vignettes

- `getting-started`:

  Quick start with no dependencies beyond base R.

- `stan-workflow`:

  Full factor analysis example in Stan with sign switching
  demonstration.

## References

Xiao, X. and Rabe-Hesketh, S. (2026). A Parameterization-Invariant DIC.
*arXiv preprint* arXiv:2605.27844.

## See also

Useful links:

- <https://github.com/DoriaXiao/DIC_i>

- Report bugs at <https://github.com/DoriaXiao/DIC_i/issues>

## Author

**Maintainer**: Xingyao Xiao <xiaoxg@berkeley.edu>

Authors:

- Sophia Rabe-Hesketh <sophiarh@berkeley.edu>
