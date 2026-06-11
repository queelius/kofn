# Fit k-out-of-n system MLE under Scheme 1 (periodic inspection)

Returns a closure that fits the model to Scheme 1 data using multi-start
optimization with L-BFGS-B and Nelder-Mead fallback.

## Usage

``` r
fit_scheme1(model, ...)
```

## Arguments

- model:

  A `kofn` model object (exponential or Weibull).

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par0 = NULL, n_starts = 5L)` that returns a
`fisher_mle` object (from the likelihood.model package).

## Details

The solver uses L-BFGS-B as the primary optimization method with
positivity constraints, falling back to Nelder-Mead on the log-parameter
scale if L-BFGS-B fails to converge.

Standard errors are computed from the numerical Hessian at the MLE.

## Examples

``` r
# \donttest{
model <- kofn(k = 2, m = 2, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
set.seed(42)
df <- rdata_scheme1(model)(c(1, 2), n = 50, delta = 1.0)
#> Error: object 'model' not found
result <- fit_scheme1(model)(df)
#> Error: object 'model' not found
coef(result)
#> Error: object 'result' not found
# }
```
