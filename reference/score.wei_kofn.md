# Score function for Weibull k-out-of-n system

Returns a closure that computes the score (gradient of log-likelihood)
via numerical differentiation using
[`numDeriv::grad`](https://rdrr.io/pkg/numDeriv/man/grad.html).

## Usage

``` r
# S3 method for class 'wei_kofn'
score(model, ...)
```

## Arguments

- model:

  A `wei_kofn` object.

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par)` returning a numeric vector (the gradient
of the log-likelihood).

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_weibull())
#> Error in dfr_weibull(): could not find function "dfr_weibull"
sc <- score(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata(model)(c(1.5, 2.0, 2.0, 3.0), n = 30)
#> Error: object 'model' not found
sc(df, c(1.5, 2.0, 2.0, 3.0))
#> Error in sc(df, c(1.5, 2, 2, 3)): could not find function "sc"
```
