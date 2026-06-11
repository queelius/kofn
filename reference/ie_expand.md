# Inclusion-exclusion expansion of a product of CDFs

For a set of exponential rates `lam`, computes the inclusion-exclusion
expansion of \\\prod\_{i} (1 - e^{-\lambda_i t})\\.

## Usage

``` r
ie_expand(lam)
```

## Arguments

- lam:

  Numeric vector of positive rates (one per component).

## Value

A list with elements:

- sign:

  Integer vector of signs (\\+1\\ or \\-1\\).

- rate_sum:

  Numeric vector of summed rates for each term.

## Details

Returns sign and rate-sum vectors such that: \$\$\prod_i (1 -
e^{-\lambda_i t}) = \sum_k \text{sign}\_k \cdot e^{-\text{rate\\sum}\_k
\cdot t}\$\$

The number of terms is \\2^m\\ where \\m\\ is `length(lam)`.

## Examples

``` r
ie <- ie_expand(c(1, 2))
# Product: (1 - exp(-t)) * (1 - exp(-2t))
# = 1 - exp(-t) - exp(-2t) + exp(-3t)
ie$sign      # c(1, -1, -1, 1)
#> [1]  1 -1 -1  1
ie$rate_sum  # c(0, 1, 2, 3)
#> [1] 0 1 2 3
```
