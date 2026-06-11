# Closed-form integral of w_j(t) over an interval

Computes \\\int_a^b w_j(t)\\ dt\\ in closed form using the
inclusion-exclusion expansion. Each term integrates as: \$\$\int_a^b
e^{-r t}\\ dt = \frac{e^{-ra} - e^{-rb}}{r}\$\$

## Usage

``` r
w_j_integral(a, b, par, j)
```

## Arguments

- a:

  Lower bound of integration (non-negative numeric scalar).

- b:

  Upper bound of integration (positive numeric scalar, `b > a`).

- par:

  Numeric vector of rates (length m).

- j:

  Component index (integer, 1-based).

## Value

Scalar integral value.

## See also

[`w_j_exact()`](https://queelius.github.io/kofn/reference/w_j_exact.md)
for the pointwise evaluation.

## Examples

``` r
# Integral of w_1(t) from 0 to 1, rates = c(1, 2)
w_j_integral(0, 1, c(1, 2), j = 1)
#> [1] 0.3153829
```
