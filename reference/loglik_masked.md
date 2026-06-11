# Masked k-out-of-n log-likelihood

Computes the log-likelihood for a k-out-of-n system with candidate
failed sets (generalized masking). The data frame must contain columns
for system lifetime, observation type, and Boolean candidate set
indicators (`c1, c2, ..., cm`).

## Usage

``` r
loglik_masked(model, candset = "c", ...)
```

## Arguments

- model:

  A `kofn` model object.

- candset:

  Prefix for candidate set columns (default `"c"`).

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par)` returning a scalar log-likelihood.

## Details

At system failure, k components have failed. The candidate set `C_i` is
a superset of the true failed set. The likelihood marginalizes over all
\\\binom{\|C_i\|}{k}\\ possible k-subsets.

For series systems (`k = 1`), this reduces to the standard masked
cause-of-failure likelihood \\\sum\_{j \in C_i} w_j(t_i)\\. For parallel
systems (`k = m`), the candidate set must be `{1, ..., m}` and the
likelihood equals the system density.

The function uses the exponential or Weibull distribution functions
directly (not the IE expansion), so it works for any family supported by
[`parse_params`](https://queelius.github.io/kofn/reference/parse_params.md).

## Examples

``` r
# 2-out-of-4 system with candidate failed sets
model <- kofn(k = 2, m = 4, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
ll <- loglik_masked(model)
#> Error: object 'model' not found

# Manually construct masked data: k=2 failed, C = {1,2,3}
df <- data.frame(
  t = 1.5, omega = "exact",
  c1 = TRUE, c2 = TRUE, c3 = TRUE, c4 = FALSE
)
ll(df, c(1, 0.8, 0.6, 0.4))
#> Error in ll(df, c(1, 0.8, 0.6, 0.4)): could not find function "ll"
```
