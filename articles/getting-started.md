# Getting started with kofn

``` r
library(kofn)
library(flexhaz)
```

## What `kofn` does

`kofn` is a maximum likelihood estimation toolkit for **k-out-of-n
systems**: reliability systems that fail when `k` of `m` components have
failed. `k = 1` is a series system (first failure ends it), `k = m` is a
parallel system (all must fail), and intermediate values interpolate
between them.

The package focuses exclusively on inference: log-likelihood, score,
Hessian, fit, observation schemes (right/left/interval/periodic), masked
cause-of-failure, and Fisher information comparison. The underlying
data-generating process (DGP) and topology queries are delegated to
`dist.structure`, which kofn imports.

## Your first fit

``` r
# A 2-of-3 exponential system: system fails when 2 of 3 components fail.
model <- kofn(k = 2, m = 3, component = dfr_exponential())

# Generate 200 complete (Scheme 0) observations.
gen <- rdata(model)
df <- gen(theta = c(1, 2, 3), n = 100)

# Fit via maximum likelihood.
fitter <- fit(model)
result <- fitter(df, n_starts = 1)
result$converged
#> [1] TRUE
coef(result)
#> [1] 1.684274 1.684549 1.682958
```

For exponential parallel systems (k = m), only the sum of rates is
identifiable from system-level data alone; individual rates are not. The
sum typically recovers well:

``` r
sum(coef(result))
#> [1] 5.051782
sum(c(1, 2, 3))
#> [1] 6
```

## Convention

kofn uses the **:F convention**: `k` counts component failures that
trigger system failure.

- `kofn(k = 1, m)` = series (first failure ends the system).
- `kofn(k = m, m)` = parallel (system survives until the m-th failure).
- `kofn(k, m)` for intermediate k = k-out-of-m systems.

(dist.structure internally uses the :G convention, which is the dual.
kofn handles the conversion for you; you never need to think about :G.)

## Observation schemes

Real experiments rarely observe system failure exactly. kofn provides
composable observation functors for censoring and periodic inspection:

``` r
# Right-censored: systems that fail after tau = 2 get recorded at t = 2
df_cens <- gen(c(1, 2, 3), n = 100,
               observe = observe_right_censor(tau = 2))
table(df_cens$omega)
#> 
#> exact 
#>   100
```

The log-likelihood handles mixed observation types automatically:

``` r
ll <- loglik(model)
ll(df_cens, c(1, 2, 3))  # finite; right-censored obs contribute log S(t)
#> [1] -26.65438
```

Other functors: `observe_left_censor`, `observe_interval_censor`,
`observe_periodic`, `observe_mixture`. See the observation-schemes
vignette.

## Weibull components

Weibull component lifetimes work analogously, with parameters
interleaved as `(shape_1, scale_1, shape_2, scale_2, ...)`:

``` r
model_wei <- kofn(k = 2, m = 2, component = dfr_weibull(), method = "em")
gen_wei <- rdata(model_wei)
df_wei <- gen_wei(theta = c(1.5, 2.0, 2.0, 3.0), n = 100)
result_wei <- fit(model_wei)(df_wei, n_starts = 1)
result_wei$shapes
#> [1] 1.892926 2.470886
result_wei$scales
#> [1] 1.485970 3.640238
```

For parallel Weibull systems (`k = m`), `method = "em"` uses an EM
algorithm treating the identity of the last-failing component as latent.
For general k, use `method = "mle"` (direct L-BFGS-B).

## Periodic inspection (Scheme 1)

When only periodic snapshots are available, use
`loglik_scheme1`/`fit_scheme1`:

``` r
s1gen <- rdata_scheme1(model)
df_s1 <- s1gen(theta = c(1, 2, 3), n = 200, delta = 0.5)
ll_s1 <- loglik_scheme1(model)
ll_s1(df_s1, c(1, 2, 3))
#> [1] -731.7861
```

## Masked cause-of-failure

When you observe system failure and a candidate set (superset of the
truly-failed components) rather than knowing which components failed:

``` r
rgen <- rdata_masked(model)
df_masked <- rgen(theta = c(1, 2, 3), n = 100, p_mask = 0.3)
ll_masked <- loglik_masked(model)
ll_masked(df_masked, c(1, 2, 3))
#> [1] -58.58316
```

## Fisher information comparison

How much precision do you gain from each observation scheme? kofn
provides a built-in comparison:

``` r
res <- compare_fisher_info(
  rates = c(0.5, 0.3), n = 50, delta = 1.0, n_rep = 10,
  component = dfr_exponential()
)
res$median_det
```

This runs `n_rep` replications of each scheme and reports the median
determinant of the observed information matrix per scheme.

## What’s next

See the other vignettes for deeper coverage:

- **[Exponential
  parallel](https://queelius.github.io/kofn/articles/exponential-parallel.md)**:
  the IE expansion fast path for parallel exponentials.
- **[Observation
  schemes](https://queelius.github.io/kofn/articles/observation-schemes.md)**:
  composable observation functors for complex monitoring setups.
- **[Weibull
  EM](https://queelius.github.io/kofn/articles/weibull-em.md)**: the EM
  algorithm details.
- **[Periodic
  inspection](https://queelius.github.io/kofn/articles/periodic-inspection.md)**:
  Scheme 1 workflows.
- **[dist.structure
  integration](https://queelius.github.io/kofn/articles/dist-structure-integration.md)**:
  how kofn delegates DGP and topology to dist.structure; migration notes
  for v0.2.0 users.
- **[Ecosystem](https://queelius.github.io/kofn/articles/ecosystem.md)**:
  how kofn fits with `likelihood.model`, `algebraic.mle`, and
  `dist.structure`.
