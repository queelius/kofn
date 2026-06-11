# System survival function for Weibull parallel system (internal)

Computes S_sys(t) = 1 - prod_j F_j(t) for a parallel system.

## Usage

``` r
weibull_S_sys(t, shapes, scales)
```

## Arguments

- t:

  Numeric scalar. Time point.

- shapes:

  Numeric vector. Weibull shape parameters.

- scales:

  Numeric vector. Weibull scale parameters.

## Value

Numeric scalar S_sys(t).
