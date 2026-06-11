# Fit Weibull k-out-of-n system model

Returns a closure that fits the model to data. Two methods are
available:

## Usage

``` r
# S3 method for class 'wei_kofn'
fit(object, ...)
```

## Arguments

- object:

  A `wei_kofn` object.

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par0, ...)` that returns a `fisher_mle` object
(from the likelihood.model package).

## Details

- `"em"`:

  EM algorithm for parallel systems (k = m only). Treats the identity of
  the last-failing component as latent data. Uses truncated Weibull
  moments for the E-step and profile optimization over shape with
  closed-form scale for the M-step.

- `"mle"`:

  Direct MLE via L-BFGS-B with Nelder-Mead fallback. Works for any
  k-out-of-n structure.

The returned solver accepts these additional arguments:

- `par0`:

  Initial parameter vector (length 2\*m, interleaved shape/scale). If
  `NULL`, method-of-moments initialization is used.

- `n_starts`:

  Number of multi-start attempts (default: 5).

- `tol`:

  Convergence tolerance for EM (default: 1e-8).

- `maxiter`:

  Maximum EM iterations (default: 2000).

- `shape_bounds`:

  Bounds for shape optimization in EM (default: `c(0.05, 15)`).

- `verbose`:

  Print iteration progress (default: FALSE).

## Examples

``` r
# \donttest{
model <- kofn(k = 2, m = 2, component = dfr_weibull())
#> Error in dfr_weibull(): could not find function "dfr_weibull"
set.seed(42)
df <- rdata(model)(c(1.5, 2.0, 2.0, 3.0), n = 50)
#> Error: object 'model' not found
result <- fit(model)(df)
#> Error: object 'model' not found
coef(result)
#> Error: object 'result' not found
# }
```
