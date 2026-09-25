# Compare DIC_i Across Models

Compares DIC_i values from multiple fitted models, returning a summary
table sorted by DIC_i (lowest = preferred).

## Usage

``` r
compare_dic_i(...)
```

## Arguments

- ...:

  Named `DICi` objects (output of
  [`compute_dic_i()`](https://doriaxiao.github.io/DIC_i/reference/compute_dic_i.md)).

## Value

A data frame with columns `model`, `dic_i`, `p_v`, `e_d`, `delta_dic_i`
(difference from the best model).

## Examples

``` r
set.seed(42)
ll_1 <- matrix(rnorm(500 * 50, -2.0, 0.5), 500, 50)
ll_2 <- matrix(rnorm(500 * 50, -2.1, 0.5), 500, 50)

fit1 <- compute_dic_i(ll_1)
fit2 <- compute_dic_i(ll_2)

compare_dic_i(model_1 = fit1, model_2 = fit2)
#>     model    dic_i      p_v      e_d delta_dic_i
#> 1 model_1 224.1131 24.03493 200.0782      0.0000
#> 2 model_2 235.7243 25.52613 210.1982     11.6112
```
