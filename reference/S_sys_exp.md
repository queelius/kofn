# System survival function for exponential parallel systems

Computes \\S\_{sys}(t) = 1 - F\_{sys}(t)\\ for a parallel system with
exponential components.

## Usage

``` r
S_sys_exp(t, par)
```

## Arguments

- t:

  Scalar time point (non-negative numeric).

- par:

  Numeric vector of rates (length m).

## Value

Scalar survival probability \\P(T\_{sys} \> t)\\.

## See also

[`F_sys_exp()`](https://queelius.github.io/kofn/reference/F_sys_exp.md)
for the CDF.

## Examples

``` r
S_sys_exp(1, c(1, 2))  # P(T_sys > 1) for 2-component parallel
#> [1] 0.4534277
```
