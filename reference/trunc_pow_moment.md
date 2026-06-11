# Scalar truncated power moment of the Weibull distribution

Convenience wrapper around
[`trunc_pow_moment_vec()`](https://queelius.github.io/kofn/reference/trunc_pow_moment_vec.md)
for a single truncation point. Computes \\E\[T^k \| T \< t\]\\ for \\T
\sim \text{Weibull}(\alpha, \beta)\\.

## Usage

``` r
trunc_pow_moment(k, t, alpha, beta)
```

## Arguments

- k:

  Numeric scalar. Power parameter (can be non-integer).

- t:

  Numeric scalar. Truncation point (positive).

- alpha:

  Numeric scalar. Weibull shape parameter.

- beta:

  Numeric scalar. Weibull scale parameter.

## Value

Numeric scalar \\E\[T^k \| T \< t\]\\.

## See also

[`trunc_pow_moment_vec()`](https://queelius.github.io/kofn/reference/trunc_pow_moment_vec.md)
for the vectorized version.

## Examples

``` r
# E[T^1 | T < 2] for Weibull(shape=1.5, scale=1)
trunc_pow_moment(k = 1, t = 2, alpha = 1.5, beta = 1)
#> [1] 0.8067229
```
