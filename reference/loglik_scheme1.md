# Log-likelihood for Scheme 1 (periodic inspection) parallel system

Returns a closure that computes the log-likelihood under Scheme 1
observation. The likelihood combines the system density at the exact
system failure time with interval-censored component contributions.

## Usage

``` r
loglik_scheme1(model, ...)
```

## Arguments

- model:

  A `kofn` model object (exponential or Weibull).

- ...:

  Additional arguments (currently unused).

## Value

A function `function(df, par)` returning a scalar log-likelihood.

## Details

The log-likelihood for observation i is: \$\$\log L_i(\theta) = \log
f\_{sys}(t_i) + \sum_j \log\[F_j(a\_{ij}^+) - F_j(a\_{ij}^-)\]\$\$

where \\f\_{sys}\\ is the parallel system density and \\\[a\_{ij}^-,
a\_{ij}^+)\\ is the inspection interval containing component j's failure
time.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
ll <- loglik_scheme1(model)
#> Error: object 'model' not found
set.seed(1)
df <- rdata_scheme1(model)(c(1, 2), n = 30, delta = 1.0)
#> Error: object 'model' not found
ll(df, c(1, 2))
#> Error in ll(df, c(1, 2)): could not find function "ll"
```
