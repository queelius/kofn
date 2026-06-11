# Hessian of the log-likelihood for exponential k-out-of-n model

Returns a closure `function(df, par)` that computes the Hessian matrix
of the log-likelihood with respect to the rate parameters.

## Usage

``` r
# S3 method for class 'exp_kofn'
hess_loglik(model, ...)
```

## Arguments

- model:

  An `exp_kofn` object created by
  [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md).

- ...:

  Additional arguments (ignored).

## Value

A function `function(df, par)` returning an `m x m` Hessian matrix.

## Details

Uses numerical differentiation via
[`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
applied to the log-likelihood closure.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
H <- hess_loglik(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata(model)(c(1, 2), n = 30)
#> Error: object 'model' not found
H(df, c(1, 2))  # 2x2 Hessian matrix
#> Error in H(df, c(1, 2)): could not find function "H"
```
