# MLE fitting for exponential k-out-of-n model

Returns a closure `function(df, par0, n_starts)` that computes the
maximum likelihood estimate of the component rate parameters.

## Usage

``` r
# S3 method for class 'exp_kofn'
fit(object, ...)
```

## Arguments

- object:

  An `exp_kofn` object created by
  [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md).

- ...:

  Additional arguments (ignored).

## Value

A function `function(df, par0 = NULL, n_starts = 5L)` returning a
fisher_mle object (from likelihood.model).

## Details

Uses multi-start optimization with L-BFGS-B as primary solver and
Nelder-Mead on the log-scale as fallback. Standard errors are computed
from the observed Fisher information (negative Hessian) at the MLE.

## Examples

``` r
# \donttest{
model <- kofn(k = 2, m = 2, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
set.seed(42)
df <- rdata(model)(c(1, 2), n = 50)
#> Error: object 'model' not found
result <- fit(model)(df)
#> Error: object 'model' not found
coef(result)
#> Error: object 'result' not found
# }
```
