# Log-likelihood for exponential k-out-of-n model

Returns a closure `function(df, par)` that computes the log-likelihood
of exponential component rates given system-level data.

## Usage

``` r
# S3 method for class 'exp_kofn'
loglik(model, ...)
```

## Arguments

- model:

  An `exp_kofn` object created by
  [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md).

- ...:

  Additional arguments (ignored).

## Value

A function `function(df, par)` where `df` is a data frame with columns
for lifetime, observation type, and (for interval-censored) upper
bounds, and `par` is a numeric vector of `m` component rates.

## Details

For parallel systems (`k = m`), the closed-form IE expansion is used,
supporting all four observation types: exact, right, left, and interval
censored. For general k-out-of-n systems, the per-observation
contributions are computed via `dist.structure::exp_kofn(k, par)` and
the algebraic.dist generics `density`, `surv`, and `cdf`.

## Examples

``` r
model <- kofn(k = 3, m = 3, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
ll <- loglik(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata(model)(c(1, 2, 3), n = 30)
#> Error: object 'model' not found
ll(df, c(1, 2, 3))
#> Error in ll(df, c(1, 2, 3)): could not find function "ll"
```
