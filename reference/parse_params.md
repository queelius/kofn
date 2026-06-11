# Parse a flat parameter vector into shapes and scales

Unified extraction of Weibull shape/scale parameters from the flat
vectors used by optimizers. For exponential components, shapes are all 1
and scales are 1/rate.

## Usage

``` r
parse_params(par, m, component)
```

## Arguments

- par:

  Numeric parameter vector. Length `m` for exponential, `2*m` for
  Weibull.

- m:

  Number of components.

- component:

  A `dfr_dist` prototype (e.g.
  [`dfr_exponential()`](https://queelius.github.io/flexhaz/reference/dfr_exponential.html)
  or
  [`dfr_weibull()`](https://queelius.github.io/flexhaz/reference/dfr_weibull.html))
  that determines how `par` is laid out.

## Value

A list with components `shapes` and `scales`.
