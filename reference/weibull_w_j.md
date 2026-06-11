# Component weight for Weibull parallel system (internal)

Computes w_j(t) = f_j(t) \* prod_i != j F_i(t), the contribution of
component j to the parallel system density at time t.

## Usage

``` r
weibull_w_j(t, shapes, scales, j)
```

## Arguments

- t:

  Numeric scalar. Time point (positive).

- shapes:

  Numeric vector. Weibull shape parameters (length m).

- scales:

  Numeric vector. Weibull scale parameters (length m).

- j:

  Integer scalar. Component index (1 to m).

## Value

Numeric scalar w_j(t).
