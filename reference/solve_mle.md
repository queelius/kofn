# Default MLE solver for positive parameters

Multi-start optimization via compositional.mle: L-BFGS-B with positivity
constraints as the primary solver, Nelder-Mead on the log-parameter
scale as fallback. Runs from `n_starts` random perturbations of `par0`
and returns the best result.

## Usage

``` r
solve_mle(neg_ll, par0, n_par, n_starts = 5L, nobs = NULL)
```

## Arguments

- neg_ll:

  Function of `par` only: guarded negative log-likelihood (returns
  `.Machine$double.xmax / 2` on failure).

- par0:

  Numeric vector of initial parameter values.

- n_par:

  Integer number of parameters.

- n_starts:

  Integer number of random restarts (default 5).

- nobs:

  Integer number of observations.

## Value

An `mle_numerical` result object (from algebraic.mle) with
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html), etc. Returns `NULL`
if all starts fail.

## Examples

``` r
# Fit a simple exponential rate from positive data
x <- rexp(100, rate = 0.5)
neg_ll <- function(rate) -sum(dexp(x, rate, log = TRUE))
fit <- solve_mle(neg_ll, par0 = 1, n_par = 1, nobs = length(x))
#> Warning: NaNs produced
#> Warning: one-dimensional optimization by Nelder-Mead is unreliable:
#> use "Brent" or optimize() directly
#> Warning: NaNs produced
#> Warning: one-dimensional optimization by Nelder-Mead is unreliable:
#> use "Brent" or optimize() directly
#> Warning: NaNs produced
#> Warning: one-dimensional optimization by Nelder-Mead is unreliable:
#> use "Brent" or optimize() directly
coef(fit)
#> [1] 0.4828732
```
