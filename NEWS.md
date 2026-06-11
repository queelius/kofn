# kofn 0.4.0

## Breaking changes

The `kofn()` constructor now takes a `component` argument (a `dfr_dist`
prototype from `flexhaz`) instead of a stringy `family` argument:

```r
# Before (0.3.x):
kofn(k = 2, m = 3, family = "exponential")
kofn(k = 2, m = 3, family = "weibull")

# After (0.4.0):
kofn(k = 2, m = 3, component = dfr_exponential())
kofn(k = 2, m = 3, component = dfr_weibull())
```

Rationale: object-based component specification composes with the rest
of the ecosystem (`flexhaz`, `serieshaz`, `maskedhaz`) and makes the
supported-family table a property of class dispatch rather than a
parallel lookup in `match.arg`. Adding future supported component
families is now a one-line change.

The same rename applies to `compare_fisher_info()`:

```r
# Before:
compare_fisher_info(rates = c(1, 2), family = "exponential", ...)

# After:
compare_fisher_info(rates = c(1, 2), component = dfr_exponential(), ...)
```

Internal helpers `parse_params()`, `kofn_dgp()`, and `kofn_components()`
take `component` rather than `family`. These were `@keywords internal`
and not part of the public API, but downstream code that touched them
via `:::` will need the same rename.

## New internal helpers

* `n_par_kofn(m, component)`: expected parameter count for a homogeneous
  kofn system, replacing inline `switch(family, ...)` logic.
* `kofn_subclass(component)`: maps a `dfr_dist` prototype to the
  corresponding kofn subclass (`exp_kofn` or `wei_kofn`).

## Dependencies

* Added `flexhaz` (>= 0.5.2) to Imports. The `dfr_dist` prototype
  constructors provide family-specific S3 subclasses (`dfr_exponential`,
  `dfr_weibull`) that `kofn_subclass()` dispatches on.

## Bug fixes

* `compare_fisher_info()`: corrected the Weibull complete-data
  (Scheme 2) Fisher information formula. The shape-shape entry had a
  spurious leading 1 (about 55 percent too large) and the shape-scale
  cross term was divided by `alpha * beta` instead of `beta`. Both
  errors vanish at shape 1 (the exponential boundary), which is why
  exponential results were unaffected. The analytic block now lives in
  an internal helper verified against numerically integrated expected
  information in the test suite. Scheme 0/1 results (numerical
  Hessians) were never affected; for Weibull components, Scheme 2
  determinants and the efficiency ratios involving Scheme 2 change.
* Interval-censored rows in a data frame lacking the upper-bound
  column (`t_upper` by default) now fail with an informative error
  in all likelihood paths. Previously the exponential parallel fast
  path crashed with "missing value where TRUE/FALSE needed" and the
  general-k path silently returned `-Inf`.

# kofn 0.3.1

## Documentation

* Vignettes restructured around the post-refactor package scope.
  Deleted outdated vignettes (`censoring-framework`, `kofn-systems`,
  `general-kofn`, `symmetry-breaking`) that documented removed
  topology / general-system APIs; added:
  - `getting-started`: quick tour of the kofn model, fit, and
    observation schemes.
  - `dist-structure-integration`: how kofn delegates DGP and topology
    to dist.structure, convention conversion, migration notes for
    v0.2.0 users.
* Added `_pkgdown.yml` configuring navbar, article menus, and a
  reference index grouped by theme.
* Retained and verified: `exponential-parallel`, `observation-schemes`,
  `weibull-em`, `periodic-inspection`, `ecosystem`.

# kofn 0.3.0

## Major refactor: dist.structure adoption

kofn now delegates the data-generating process (DGP), topology, and
system-level distribution generics to `dist.structure`. The package
focuses exclusively on inference for k-out-of-n systems: log-likelihood,
score, Hessian, fit, observation schemes (right/left/interval/periodic),
masked cause-of-failure, periodic-inspection (Scheme 1), and Fisher
information comparison.

### Removed

- `coherent_system()`, `parallel_system()`, `series_system()`,
  `kofn_system()`, `bridge_system()`, `consecutive_k_system()`: use
  the corresponding constructors in `dist.structure` (`coherent_dist`,
  `parallel_dist`, `series_dist`, `kofn_dist`, `bridge_dist`,
  `consecutive_k_dist`).
- `system_signature()`, `system_lifetime()`, `system_censoring()`,
  `critical_states()`, `min_paths()`, `min_cuts()`, `phi()`,
  `min_cuts_from_paths()`, `minimize_sets()`: use `dist.structure`
  equivalents (registered on every `dist_structure` subclass).
- `loglik_system()`, `fit_system()`, `rdata_system()`: replaced by
  the `loglik`/`fit`/`rdata` S3 methods on `kofn` model objects (for
  k-of-n) plus the dist generics on `dist.structure::exp_kofn`,
  `wei_kofn`, `coherent_dist` (for general topologies).
- `f_sys_general()`, `S_sys_general()`: use
  `algebraic.dist::density()` and `surv()` on a `dist.structure`
  object.
- `make_dists()`, `exp_dist()`, `weibull_dist()` (the lightweight
  `$pdf`/`$cdf`/`$surv` list interface): replaced by
  `algebraic.dist::exponential()` and `weibull_dist()` plus the
  algebraic.dist generics.

### Changed

- `kofn()` constructor no longer accepts a `system = ...` argument.
  kofn is now exclusively for k-of-n estimation. Users with non-k-of-n
  topologies should use `dist.structure::coherent_dist()` directly.
- `loglik.exp_kofn` and `loglik.wei_kofn` general-k path now delegates
  per-observation contributions to `dist.structure::exp_kofn(k_dist,
  par)` and `wei_kofn(k_dist, shapes, scales)` via the algebraic.dist
  `density`, `surv`, and `cdf` generics. Parallel-fast-path (k=m) IE
  expansion is preserved as the closed-form path.
- `loglik_scheme1` system density now uses
  `density(dist.structure::exp_kofn(...))` / `wei_kofn(...)` instead
  of the local `f_sys_general` engine.
- `rdata.exp_kofn`, `rdata.wei_kofn`, `rdata_masked`,
  `rdata_scheme1` use a kofn-internal specialized `kofn_censoring(k,
  times)` helper (the k-th order statistic) instead of the general
  `system_censoring` algorithm.
- `assumptions.exp_kofn` and `assumptions.wei_kofn` no longer mention
  "general coherent system" since kofn is k-of-n only.

### Internal

- New file `R/internal_topology.R` with private helpers
  `kofn_censoring()`, `kofn_dgp()`, `kofn_components()`. Uses the
  k_dist = m - k_kofn + 1 conversion between kofn's :F convention
  (k_kofn = number of failures triggering system failure) and
  dist.structure's :G convention (k_dist = number of functioning
  components required).
- Deleted `R/coherent_system.R` (599 lines), `R/system_density.R`
  (336 lines), and `R/dist.R` (118 lines).
- Package shrank from 4312 lines to roughly 2750 lines (~36%
  reduction) in the R/ directory.

### Added

- `Imports: dist.structure, algebraic.dist`.
- `ncomponents` re-exported from dist.structure (was a local generic).

# kofn 0.2.0

* Generalized masking: `loglik_masked()` and `rdata_masked()` extend masked
  cause-of-failure from series (k=1) to arbitrary k-out-of-m systems.

* Periodic inspection (Scheme 1): `loglik_scheme1()` and `fit_scheme1()`
  support general k via the critical-state density engine.

* Performance: eagerly precomputed critical/functioning states in
  `coherent_system()` and vectorized `f_sys_general()` give ~70x speedup
  for Scheme 1 likelihood evaluation.

* Migrated optimization to the `compositional.mle` package. The internal
  helper `multistart_mle()` was removed and replaced by a private
  `solve_mle()` wrapper around `compositional.mle::lbfgsb()` with a
  Nelder-Mead fallback in log-parameterization.

* New vignettes: symmetry-breaking comparison, periodic inspection,
  observation schemes, general k-out-of-m.

# kofn 0.1.0

* Initial release with exponential and Weibull parallel system MLE.
* EM algorithm for Weibull parallel systems.
* Coherent system infrastructure via minimal path/cut sets.
* Fisher information comparison across observation schemes.
