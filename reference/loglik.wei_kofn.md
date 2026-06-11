# Log-likelihood for Weibull k-out-of-n system

Returns a closure `function(df, par)` that computes the log-likelihood
for a Weibull k-out-of-n system given data and parameters.

## Usage

``` r
# S3 method for class 'wei_kofn'
loglik(model, ...)
```

## Arguments

- model:

  A `wei_kofn` object.

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par)` returning a scalar log-likelihood.

## Details

For parallel systems (k = m), uses the direct Weibull system density
\\f\_{sys}(t) = \sum_j f_j(t) \prod\_{i \ne j} F_i(t)\\. Supports
`"exact"` and `"right"` observation types. Left and interval censoring
are not supported for Weibull (no IE expansion); use the exponential
model or Scheme 1 for those cases.

For general k, delegates per-observation to
`dist.structure::wei_kofn(k, shapes, scales)` via the algebraic.dist
generics `density`, `surv`, and `cdf`.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_weibull())
#> Error in dfr_weibull(): could not find function "dfr_weibull"
ll <- loglik(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata(model)(c(1.5, 2.0, 2.0, 3.0), n = 30)
#> Error: object 'model' not found
ll(df, c(1.5, 2.0, 2.0, 3.0))
#> Error in ll(df, c(1.5, 2, 2, 3)): could not find function "ll"
```
