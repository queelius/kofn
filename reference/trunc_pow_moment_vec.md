# Vectorized truncated power moment of the Weibull distribution

Computes \\E\[T^k \| T \< t_i\]\\ for \\T \sim \text{Weibull}(\alpha,
\beta)\\ simultaneously for a vector of truncation points.

## Usage

``` r
trunc_pow_moment_vec(k, v_max_vec, alpha, beta)
```

## Arguments

- k:

  Numeric scalar. Power parameter (can be non-integer).

- v_max_vec:

  Numeric vector. Precomputed values of \\(t_i/\beta)^\alpha\\.

- alpha:

  Numeric scalar. Weibull shape parameter.

- beta:

  Numeric scalar. Weibull scale parameter.

## Value

Numeric vector of \\E\[T^k \| T \< t_i\]\\ values, same length as
`v_max_vec`.

## Details

Uses the substitution \\U = (T/\beta)^\alpha \sim \text{Exp}(1)\\, so
that \\T^k = \beta^k U^{k/\alpha}\\ and the truncation \\T \< t\\
becomes \\U \< v\_{\max}\\. The result is: \$\$E\[T^k \| T \< t\] =
\beta^k \frac{\gamma(k/\alpha + 1, v\_{\max})}{1 - e^{-v\_{\max}}}\$\$
where \\\gamma(s, x)\\ is the lower incomplete gamma function.

Computation is done in log-space for numerical stability.

For very small \\v\_{\max}\\ (below `1e-12`), the asymptotic
approximation \\E\[T^k \| T \< t\] \approx t^k / (k/\alpha + 1)\\ is
used to avoid numerical issues.

## Examples

``` r
# E[T^2 | T < t] for Weibull(shape=2, scale=1) at t = 0.5, 1.0, 2.0
v_max <- (c(0.5, 1.0, 2.0) / 1)^2
trunc_pow_moment_vec(k = 2, v_max_vec = v_max, alpha = 2, beta = 1)
#> [1] 0.1197971 0.4180233 0.9253706
```
