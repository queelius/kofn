# Compare Fisher information across observation schemes

For a given parameter configuration, computes the observed Fisher
information under Scheme 0 (system only), Scheme 1 (periodic
inspection), and Scheme 2 (complete component data) via simulation.

## Usage

``` r
compare_fisher_info(
  shapes = NULL,
  scales = NULL,
  rates = NULL,
  n = 200L,
  delta = 1,
  n_rep = 50L,
  component = dfr_exponential()
)
```

## Arguments

- shapes:

  Numeric vector. Weibull shape parameters (Weibull only, or `NULL` for
  exponential).

- scales:

  Numeric vector. Weibull scale parameters (Weibull only, or `NULL` for
  exponential).

- rates:

  Numeric vector. Exponential rate parameters (exponential only,
  alternative to `shapes`/`scales`).

- n:

  Integer. Sample size per replicate.

- delta:

  Numeric scalar. Inspection interval width for Scheme 1.

- n_rep:

  Integer. Number of Monte Carlo replicates.

- component:

  A `dfr_dist` prototype (e.g.
  [`dfr_exponential()`](https://queelius.github.io/flexhaz/reference/dfr_exponential.html)
  or
  [`dfr_weibull()`](https://queelius.github.io/flexhaz/reference/dfr_weibull.html))
  specifying the component family.

## Value

A list with components:

- `scheme0_det`:

  Numeric vector of Scheme 0 information determinants (length `n_rep`).

- `scheme1_det`:

  Numeric vector of Scheme 1 information determinants.

- `scheme2_det`:

  Numeric vector of Scheme 2 (complete data) information determinants.

- `median_det`:

  Named numeric vector of median determinants across schemes.

- `efficiency_01`:

  Ratio of median Scheme 0 to Scheme 1 determinant (\< 1 means Scheme 1
  is more informative).

- `efficiency_02`:

  Ratio of median Scheme 0 to Scheme 2 determinant.

- `efficiency_12`:

  Ratio of median Scheme 1 to Scheme 2 determinant.

- `n_valid`:

  Named integer vector of valid (non-NA) replicate counts per scheme.

## Details

The determinant of the observed information matrix is used as a scalar
summary. Efficiency ratios indicate relative information content: a
ratio less than 1 means the denominator scheme carries more information.

For Scheme 2 (complete data), the Fisher information is computed
analytically:

- Exponential: \\I\_{jj} = n / \lambda_j^2\\ (diagonal).

- Weibull: Block-diagonal with per-component 2x2 Fisher information
  matrices using the standard Weibull FIM formulas.

For Schemes 0 and 1, the observed information is computed numerically
via [`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
of the negative log-likelihood evaluated at the true parameter values.

## Examples

``` r
# \donttest{
set.seed(1)
result <- compare_fisher_info(rates = c(1, 2), n = 50, delta = 1.0, n_rep = 5)
result$median_det
#>   scheme0   scheme1   scheme2 
#>  120.7422 1093.5996  625.0000 
result$efficiency_01  # Scheme 0 vs Scheme 1
#> [1] 0.110408
# }
```
