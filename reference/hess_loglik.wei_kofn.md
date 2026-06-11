# Hessian of log-likelihood for Weibull k-out-of-n system

Returns a closure that computes the Hessian matrix of the log-likelihood
via numerical differentiation using
[`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html).

## Usage

``` r
# S3 method for class 'wei_kofn'
hess_loglik(model, ...)
```

## Arguments

- model:

  A `wei_kofn` object.

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par)` returning a square matrix (the Hessian of
the log-likelihood).

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_weibull())
#> Error in dfr_weibull(): could not find function "dfr_weibull"
H <- hess_loglik(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata(model)(c(1.5, 2.0, 2.0, 3.0), n = 30)
#> Error: object 'model' not found
H(df, c(1.5, 2.0, 2.0, 3.0))  # 4x4 Hessian matrix
#> Error in H(df, c(1.5, 2, 2, 3)): could not find function "H"
```
