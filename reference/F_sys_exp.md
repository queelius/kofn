# System CDF for exponential parallel systems

Computes \\F\_{sys}(t) = \prod_j (1 - e^{-\lambda_j t})\\ using the
inclusion-exclusion expansion.

## Usage

``` r
F_sys_exp(t, par)
```

## Arguments

- t:

  Scalar time point (non-negative numeric).

- par:

  Numeric vector of rates (length m).

## Value

Scalar CDF value \\P(T\_{sys} \le t)\\.

## See also

[`S_sys_exp()`](https://queelius.github.io/kofn/reference/S_sys_exp.md)
for the survival function.

## Examples

``` r
F_sys_exp(1, c(1, 2))  # P(T_sys <= 1) for 2-component parallel
#> [1] 0.5465723
```
