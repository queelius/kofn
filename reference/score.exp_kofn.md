# Score function for exponential k-out-of-n model

Returns a closure `function(df, par)` that computes the gradient of the
log-likelihood with respect to the rate parameters.

## Usage

``` r
# S3 method for class 'exp_kofn'
score(model, ...)
```

## Arguments

- model:

  An `exp_kofn` object created by
  [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md).

- ...:

  Additional arguments (ignored).

## Value

A function `function(df, par)` returning a numeric vector of length `m`
(the score vector).

## Details

Uses numerical differentiation via
[`numDeriv::grad()`](https://rdrr.io/pkg/numDeriv/man/grad.html) applied
to the log-likelihood closure.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
sc <- score(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata(model)(c(1, 2), n = 30)
#> Error: object 'model' not found
sc(df, c(1, 2))  # gradient at true parameters
#> Error in sc(df, c(1, 2)): could not find function "sc"
```
