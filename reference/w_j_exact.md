# Compute w_j(t) = f_j(t) \* prod_i != j F_i(t) for exponential components

For component j with rate \\\lambda_j\\ and other rates
\\\lambda\_{-j}\\: \$\$w_j(t) = \lambda_j e^{-\lambda_j t} \prod\_{i
\neq j} (1 - e^{-\lambda_i t})\$\$

## Usage

``` r
w_j_exact(t, par, j)
```

## Arguments

- t:

  Scalar time point (positive numeric).

- par:

  Numeric vector of rates (length m), one per component.

- j:

  Component index (integer, 1-based).

## Value

Scalar value of \\w_j(t)\\.

## Details

Uses the inclusion-exclusion expansion to express this as a finite sum
of exponentials.

## See also

[`w_j_integral()`](https://queelius.github.io/kofn/reference/w_j_integral.md)
for the closed-form integral of \\w_j\\.

## Examples

``` r
# Component 1 contribution at t = 0.5, rates = c(1, 2, 3)
w_j_exact(0.5, c(1, 2, 3), j = 1)
#> [1] 0.2978523
```
