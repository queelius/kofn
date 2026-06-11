# System density for Weibull parallel system (internal)

Computes f_sys(t) = sum_j f_j(t) \* prod_i != j F_i(t).

## Usage

``` r
weibull_f_sys(t, shapes, scales)
```

## Arguments

- t:

  Numeric scalar. Time point.

- shapes:

  Numeric vector. Weibull shape parameters.

- scales:

  Numeric vector. Weibull scale parameters.

## Value

Numeric scalar f_sys(t).
