# Vectorized truncated log-moment of the Weibull distribution

Computes \\E\[\log T \| T \< t_i\]\\ for \\T \sim \text{Weibull}(\alpha,
\beta)\\ simultaneously for a vector of truncation points.

## Usage

``` r
trunc_log_moment_vec(v_max_vec, alpha, beta)
```

## Arguments

- v_max_vec:

  Numeric vector. Precomputed values of \\(t_i/\beta)^\alpha\\.

- alpha:

  Numeric scalar. Weibull shape parameter.

- beta:

  Numeric scalar. Weibull scale parameter.

## Value

Numeric vector of \\E\[\log T \| T \< t_i\]\\ values, same length as
`v_max_vec`.

## Details

Uses the identity: \$\$E\[\log T \| T \< t\] = \log \beta +
\frac{1}{\alpha} E\[\log U \| U \< v\_{\max}\]\$\$ where \\U \sim
\text{Exp}(1)\\ and \\v\_{\max} = (t/\beta)^\alpha\\.

The inner expectation \\E\[\log U \| U \< v\_{\max}\]\\ is computed via
central finite difference of the lower incomplete gamma function at \\s
= 1\\.

For very small \\v\_{\max}\\ (below `1e-10`), the asymptotic
approximation \\E\[\log T \| T \< t\] \approx \log t - 1/\alpha\\ is
used. For large \\v\_{\max}\\ where numerical issues arise, the
untruncated moment \\E\[\log T\] = \log \beta + \psi(1)/\alpha\\ is used
as fallback, where \\\psi\\ is the digamma function.

## Examples

``` r
# E[log T | T < t] for Weibull(shape=2, scale=1) at t = 0.5, 1.0, 2.0
v_max <- (c(0.5, 1.0, 2.0) / 1)^2
trunc_log_moment_vec(v_max_vec = v_max, alpha = 2, beta = 1)
#> [1] -1.2248035 -0.6301010 -0.3088497
```
