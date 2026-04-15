# kofn Rewrite Plan (Phase 1, design captured 2026-04-15)

This document captures the design for the kofn rewrite that adopts
`dist.structure` as the canonical source for k-out-of-n DGP and
topology, leaving kofn focused on its distinctive contributions:
likelihood-based inference, observation schemes, EM for Weibull,
Fisher information comparison, and Scheme 1 / masked observations.

## Mission of post-rewrite kofn

**Maximum likelihood estimation for k-out-of-n systems.** Not a general
coherent-system package, not a topology library, just estimation for
k-of-n.

## Scope decisions

### Files to delete (superseded by dist.structure)

| File | Lines | Why |
|---|---|---|
| `R/coherent_system.R` | 599 | Topology infrastructure now lives in `dist.structure::coherent_dist` and the topology shortcuts |
| `R/system_density.R` | 336 | General system density / loglik replaced by `algebraic.dist::density(dist.structure::exp_kofn(...))` etc. |

### Files to keep (kofn-specific)

| File | Lines | Role |
|---|---|---|
| `R/kofn.R` | 179 | Likelihood model constructor (simplified) |
| `R/exp_kofn.R` | 421 | Exponential MLE: loglik, score, hess_loglik, fit, rdata |
| `R/wei_kofn.R` | 666 | Weibull MLE + EM algorithm |
| `R/observe.R` | 205 | Observation functors (right-censor, periodic, mixture, ...), distinctive |
| `R/scheme1.R` | 250 | Periodic-inspection-specific likelihood |
| `R/masked.R` | 254 | Masked-cause observations (overlap with maskedcauses to evaluate) |
| `R/fisher_info.R` | 225 | Fisher information comparison across observation schemes, distinctive |
| `R/truncated_moments.R` | 173 | Weibull EM helpers |
| `R/utils.R` | 102 | Small helpers |
| `R/reexports.R` | 31 | Re-exports |

### Files to evaluate

| File | Lines | Decision |
|---|---|---|
| `R/dist.R` | 118 | Lightweight component-distribution objects with `$pdf`/`$cdf`/`$surv`. Used internally by exp_kofn / wei_kofn / masked / scheme1. Could be replaced with `dist.structure` component access, but the closure-with-fields pattern is convenient. **Keep as internal helper for now.** |
| `R/ie_expand.R` | 154 | Inclusion-exclusion expansion for parallel exponentials. Used by exp_kofn parallel-fast-path. **Keep**, it is distinctive math used by estimation. |

## Adoption steps (Phase 2, execute in dedicated session)

### Step 1, DESCRIPTION

```
 Imports:
     stats,
     numDeriv,
     likelihood.model,
     compositional.mle,
     generics,
     dist.structure,
     algebraic.dist
```

### Step 2, add new methods to dist.structure if needed

For kofn loglik/density to delegate to dist.structure cleanly:

- `density.exp_kofn(x, ...)`: closed form via critical-state formula or
  IE-based sum over components. Required for kofn's exact-observation
  loglik path to use `algebraic.dist::density(dist.structure::exp_kofn(...))(t)`
  instead of the local `f_sys_general`.
- `density.wei_kofn(x, ...)`: same for Weibull.

Without these, kofn must compute density inline. Either approach works;
adding them to dist.structure means future packages can reuse.

### Step 3, refactor exp_kofn.R

Three loglik paths in current code:
1. **Parallel (k=m) closed-form**: uses `w_j_exact`, `S_sys_exp`,
   `w_j_integral`, these are kofn-internal; keep as-is. Optionally
   verify they agree with dist.structure::exp_kofn surv/cdf as a
   correctness check.
2. **General k loglik**: currently delegates to `loglik_system` (from
   system_density.R). Replace with: compute log-density at exact
   observations + log-survival at right-censored observations using
   `algebraic.dist::density(dist.structure::exp_kofn(k, par))(t)` and
   `algebraic.dist::surv(...)`.
3. **Score and Hessian**: numerical via `numDeriv::grad` and
   `numDeriv::hessian` on the loglik closure. Unaffected.

### Step 4, refactor wei_kofn.R

Similar treatment. The EM algorithm is preserved (it is distinctive and
correct). The general-k loglik path replaces `loglik_system` calls with
`dist.structure::wei_kofn(k, shapes, scales)` based survival/density.

### Step 5, refactor masked.R and scheme1.R

`make_dists`, `S_sys_general`, `f_sys_general` calls replace with:

- `make_dists(par, "exponential")` becomes `dist.structure::exp_kofn(k, par)$components`
  (or just construct components directly via algebraic.dist)
- `S_sys_general(t, sys, dists)` becomes `algebraic.dist::surv(dist.structure::exp_kofn(k, par))(t)`
- `f_sys_general(t, sys, dists)` becomes `algebraic.dist::density(...)(t)`

If density methods are not added to dist.structure, these can be
inlined as small kofn helpers.

### Step 6, update tests

Many tests reference removed exports (parallel_system, series_system,
kofn_system, bridge_system, ie_expand, make_dists, f_sys_general,
S_sys_general, system_signature, etc.). Update to:

- `parallel_system(m)` becomes `dist.structure::parallel_dist(comps)` (with
  appropriate component constructions)
- `series_system(m)` becomes `dist.structure::series_dist(comps)`
- `kofn_system(k, m)` becomes `dist.structure::kofn_dist(k, comps)`
- `bridge_system()` becomes `dist.structure::bridge_dist(comps)`
- `system_signature(sys)` becomes `dist.structure::system_signature(sys)`

`ie_expand` is kofn-internal; tests that probe IE behavior can use
`kofn:::ie_expand` or be moved to internal-test infrastructure.

### Step 7, rewrite kofn.R constructor

Simplify: drop the optional `system = ...` argument (users who want
non-k-of-n topologies use dist.structure directly; kofn is for k-of-n).

```r
kofn <- function(k, m, family = c("exponential", "weibull"),
                 method = c("mle", "em"),
                 lifetime = "t", omega = "omega",
                 lifetime_upper = "t_upper") {
  family <- match.arg(family)
  method <- match.arg(method)
  ...
}
```

### Step 8, update NAMESPACE

Drop exports for:
- `parallel_system`, `series_system`, `kofn_system`, `bridge_system`,
  `consecutive_k_system`, `coherent_system` (use dist.structure)
- `min_cuts_from_paths`, `minimize_sets`, `min_cuts`, `critical_states`,
  `phi`, `system_signature` (use dist.structure)
- `loglik_system`, `f_sys_general`, `S_sys_general`, `make_dists`,
  `S_sys_exp`, `F_sys_exp` (internal)
- `coherent_system` print method (no longer needed)
- `exp_dist`, `solve_mle` (internal helpers)

Keep exports:
- `kofn`, `is_kofn`, `print.kofn`, `ncomponents`
- `loglik`, `score`, `hess_loglik`, `fit`, `rdata`, `assumptions`
- `loglik_masked`, `rdata_masked`, `loglik_scheme1`, `rdata_scheme1`,
  `fit_scheme1`, `fit_system`
- `compare_fisher_info`
- All `observe_*` functors

### Step 9, bump version to 0.3.0

Breaking changes warrant a minor bump. Add NEWS entry.

### Step 10, run R CMD check, fix any warnings

## Estimated effort

- Phase 2 Steps 1, 7, 8, 9: 30 min (mechanical)
- Phase 2 Steps 2, 3, 4, 5: 60-90 min (careful refactoring with test runs)
- Phase 2 Step 6: 60 min (test rewrites)
- Phase 2 Step 10: 30 min (cleanup)

Total: 3 to 4 hours of focused work. Best done in a fresh session.

## Risk areas

1. **The kofn-internal loglik for general k** depends on the
   critical-state enumeration in system_density.R. This is non-trivial
   math; replacing it with dist.structure equivalents requires either
   (a) adding density methods to dist.structure or (b) inlining a
   minimal version in kofn.

2. **The `dists` list interface** (make_dists returning list of objects
   with $pdf/$cdf/$surv) is used pervasively in masked.R and
   scheme1.R. Replacing with dist.structure components changes the
   shape; refactoring the call sites is mechanical but tedious.

3. **Test dependencies** are extensive. Updating ~15-20 tests to use
   dist.structure constructors instead of kofn's removed exports is
   the largest single chunk of work.

## Recommendation

Treat kofn rewrite as a dedicated session task. The plan above has
sufficient detail to execute in one focused sitting. Do NOT split
across sessions in flight, partial progress would leave the package
broken (failing tests, inconsistent exports).
