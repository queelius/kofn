# Mixture of observation schemes

Creates an observation functor that randomly selects from a set of
observation schemes for each observation. This models heterogeneous
monitoring environments where different units may be observed under
different protocols.

## Usage

``` r
observe_mixture(..., weights = NULL)
```

## Arguments

- ...:

  Observation functors (created by `observe_*` functions).

- weights:

  Mixing probabilities (numeric vector). If `NULL`, uniform weights are
  used. Weights are normalized to sum to 1.

## Value

Observation functor: `function(t_true)` that randomly selects a
constituent scheme and returns its output.

## Examples

``` r
obs <- observe_mixture(
  observe_right_censor(tau = 100),
  observe_periodic(delta = 10, tau = 100),
  weights = c(0.7, 0.3)
)
set.seed(42)
obs(30)
#> $t
#> [1] 30
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 40
#> 
```
