# dist.structure integration

``` r
library(kofn)
library(flexhaz)
library(dist.structure)
```

## Two packages, two missions

`kofn` and `dist.structure` have complementary roles:

| Concern                               | Package          | Generics / constructors                                                                                                                                                                                                                                      |
|---------------------------------------|------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Topology and DGP for coherent systems | `dist.structure` | `coherent_dist`, `series_dist`, `parallel_dist`, `kofn_dist`, `bridge_dist`, `exp_kofn`, `wei_kofn`, `exp_series`, `wei_series`, …; `phi`, `min_paths`, `min_cuts`, `system_signature`, `critical_states`, `reliability`, `structural_importance`, `dual`, … |
| Inference for k-out-of-n              | `kofn`           | [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md) model, `loglik`, `score`, `hess_loglik`, `fit`, `rdata`; observation schemes; EM for Weibull; Scheme 1 / masked likelihoods; Fisher info comparison                                            |

Before v0.3.0, kofn bundled both, which meant 1500+ lines of
topology/DGP code duplicated functionality that properly belongs to a
distribution-algebra package. The v0.3.0 refactor delegated that work to
dist.structure, shrinking kofn’s source by ~36% and focusing it on its
distinctive contribution: inference under varied observation schemes.

## What kofn delegates to dist.structure

The general-k log-likelihood path computes per-observation contributions
through dist.structure generics. Under the hood:

``` r
# Pseudo-code of kofn's loglik.exp_kofn general-k path:
function(df, par) {
  dgp <- dist.structure::exp_kofn(k_dist, par)  # k_dist = m - k_kofn + 1
  surv_fn <- algebraic.dist::surv(dgp)
  cdf_fn  <- algebraic.dist::cdf(dgp)
  dens_fn <- density(dgp)
  # ... then loop over observations, accumulating log contributions.
}
```

Your code calling `loglik(kofn_model)(df, par)` doesn’t see any of this;
the integration is invisible.

The parallel fast path (`k == m`) retains kofn’s inclusion-exclusion
expansion because that closed form is specific to exponential parallel
systems and faster than the general dist.structure density for that
case. Cross-validation tests confirm both paths produce identical
results.

## When you might invoke dist.structure directly

**1. Post-fit analysis on the DGP.** After you fit a kofn model, you
have the estimated parameters. Construct a `dist.structure` DGP to
exercise the full topology protocol:

``` r
model <- kofn(k = 2, m = 3, component = dfr_exponential())
df <- rdata(model)(theta = c(1, 2, 3), n = 100)
result <- fit(model)(df, n_starts = 1)
rates_hat <- coef(result)

# Build the DGP from the fitted rates.
sys_fitted <- dist.structure::exp_kofn(
  k = model$m - model$k + 1,   # convert to :G convention
  rates = rates_hat
)

# Query structural and reliability properties.
dist.structure::system_signature(sys_fitted)
#> [1] 0 1 0
dist.structure::reliability(sys_fitted, 0.8)
#> [1] 0.896
sapply(seq_len(model$m), function(j)
  dist.structure::structural_importance(sys_fitted, j))
#> [1] 0.5 0.5 0.5
```

**2. Non-k-of-n topologies.** kofn v0.3.0 removed the `system = ...`
argument from the constructor. If you need to estimate a bridge or other
arbitrary coherent system, use dist.structure directly: it provides
`coherent_dist`, `bridge_dist`, `consecutive_k_dist`, and can evaluate
loglik at user-supplied component distributions.

## Convention conversion

kofn uses the :F convention: `k_kofn` = number of failures that trigger
system failure. `k=1` is series; `k=m` is parallel.

dist.structure uses :G: `k_dist` = number of functioning components
required. `k=1` is parallel; `k=m` is series.

They convert via:

    k_dist = m - k_kofn + 1

Examples:

| kofn shape                    | kofn k              | dist.structure constructor                      |
|-------------------------------|---------------------|-------------------------------------------------|
| Series (fails on 1st failure) | `k_kofn = 1`        | `exp_kofn(k = m, ...)` or `exp_series(rates)`   |
| Parallel (fails on m-th)      | `k_kofn = m`        | `exp_kofn(k = 1, ...)` or `exp_parallel(rates)` |
| 2-of-3 (fails on 2nd of 3)    | `k_kofn = 2, m = 3` | `exp_kofn(k = 2, rates = ...)`                  |

The conversion is handled inside kofn whenever you pass a model through
loglik/fit/rdata; you only need to remember it when calling
dist.structure directly.

## Migration notes: v0.2.0 -\> v0.3.0

If you used kofn 0.2.0, here are the removed APIs and their
replacements:

| Removed (v0.2.0)                                                   | Replacement (v0.3.0)                                                                                                                                                                                                                      |
|--------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `parallel_system(m)`                                               | `dist.structure::parallel_dist(components)`                                                                                                                                                                                               |
| `series_system(m)`                                                 | `dist.structure::series_dist(components)`                                                                                                                                                                                                 |
| `kofn_system(k, m)`                                                | `dist.structure::kofn_dist(k, components)` (:G conv)                                                                                                                                                                                      |
| `bridge_system()`                                                  | `dist.structure::bridge_dist(components)`                                                                                                                                                                                                 |
| `consecutive_k_system(k, n)`                                       | `dist.structure::consecutive_k_dist(k, components)`                                                                                                                                                                                       |
| `coherent_system(min_paths, m)`                                    | `dist.structure::coherent_dist(min_paths, components, m)`                                                                                                                                                                                 |
| `kofn(system = ...)`                                               | no longer supported; use dist.structure for non-k-of-n                                                                                                                                                                                    |
| `make_dists(par, family)`                                          | `algebraic.dist::exponential(rate)` / `weibull_dist(shape, scale)`                                                                                                                                                                        |
| `loglik_system(t, sys, par, family)`                               | `-sum(log(algebraic.dist::density(dist.structure::exp_kofn(...))(t)))` (or wei_kofn, or coherent_dist)                                                                                                                                    |
| `fit_system(t, sys, family, ...)`                                  | For k-of-n: `fit(kofn(k, m, component = dfr_exponential()))(df, ...)` (or [`dfr_weibull()`](https://queelius.github.io/flexhaz/reference/dfr_weibull.html)). For general coherent: roll your own optimizer with dist.structure densities. |
| `rdata_system(sys, par, ...)`                                      | `algebraic.dist::sampler(dgp)(n)` on the appropriate dist.structure object                                                                                                                                                                |
| `system_signature(sys)`                                            | `dist.structure::system_signature(sys)`                                                                                                                                                                                                   |
| `system_censoring(sys, times)`                                     | `dist.structure::system_censoring(sys, times)` (note: returns `$system_time` and `$component_status`, not kofn’s `$T_sys`/`$critical`/`$status`)                                                                                          |
| `critical_states(sys, j)`                                          | `dist.structure::critical_states(sys, j)`                                                                                                                                                                                                 |
| `phi(sys, x)`, `min_paths`, `min_cuts`                             | `dist.structure::` equivalents                                                                                                                                                                                                            |
| `f_sys_general`, `S_sys_general`                                   | `algebraic.dist::density(dgp)` and `surv(dgp)`                                                                                                                                                                                            |
| `ie_expand`, `w_j_exact`, `w_j_integral`, `F_sys_exp`, `S_sys_exp` | Still exported from kofn (they’re kofn’s parallel-fast-path internals)                                                                                                                                                                    |

## Why this architecture

The reliability ecosystem now has a single source of truth for
coherent-system DGP and topology (dist.structure). Downstream packages
specialize:

- `serieshaz` implements the protocol on `dfr_dist_series` for arbitrary
  dynamic failure rate components.
- `maskedcauses` thin-wraps its exponential-series API over
  [`dist.structure::exp_series`](https://rdrr.io/pkg/dist.structure/man/exp_series.html).
- `kofn` provides inference for the k-out-of-n family.

Any new reliability package that needs topology or DGP machinery should
import dist.structure and add its specialized math layer on top, not
reimplement what dist.structure already provides.
