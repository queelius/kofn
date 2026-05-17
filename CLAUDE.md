# CLAUDE.md: kofn R package

## Overview

Component lifetime estimation from k-out-of-n system data, framed as a
**censoring problem**. A k-out-of-n system fails when k components have
failed; the system lifetime T is the k-th order statistic of
(T_1, ..., T_m). Observing T creates structured censoring on individual
component lifetimes.

**Key framing**: This is NOT masked cause-of-failure (that's the series
case, handled by `maskedcauses`). For parallel systems, all components
have failed by system failure time. The incomplete-data problem is
*when* each component failed, not *which* one.

## Ecosystem

- `flexhaz`: distribution prototypes (`dfr_exponential()`, `dfr_weibull()`)
- `dist.structure`: DGP and topology (kofn delegates system density here)
- `likelihood.model`: generics (loglik, score, hess_loglik, fit, rdata)
- `compositional.mle`: solver composition (%>>%, %|%, with_restarts)
- `algebraic.mle`: MLE result objects (coef, vcov, confint, AIC, BIC)
- `maskedcauses`: series systems (k=1), masked cause of failure
- `maskedhaz`: general DFR framework for series systems

## Commands

```bash
Rscript -e "devtools::load_all()"          # Interactive development
Rscript -e "devtools::document()"          # Generate roxygen2 docs
Rscript -e "devtools::test()"              # Run tests
Rscript -e "covr::package_coverage()"      # Coverage report
Rscript -e "devtools::check()"             # R CMD check
```

Note: flexhaz must be installed from source (or loaded via
`devtools::load_all`) for prototype subclasses (`dfr_exponential`,
`dfr_weibull`) to be available. The subclass dispatch is required for
`kofn_subclass()` to work.

## Architecture

### Model API: `R/kofn.R`
Single constructor `kofn(k, m, component, method)` returns an S3 object.
The `component` argument is a `dfr_dist` prototype from `flexhaz` (e.g.
`dfr_exponential()`, `dfr_weibull()`). The S3 class of the prototype
drives dispatch via `kofn_subclass()`.

Classes: `exp_kofn` (exponential) or `wei_kofn` (Weibull), inheriting
from `kofn` and `likelihood_model`.

Generic methods (from likelihood.model) return closures:
- `loglik(model)` -> `function(df, par)`
- `score(model)` -> `function(df, par)`
- `fit(model)` -> `function(df, par0, ...)`
- `rdata(model)` -> `function(theta, n, ...)`

### Topology delegation: `R/internal_topology.R`
Private helpers `kofn_censoring()`, `kofn_dgp()`, `kofn_components()`
delegate to `dist.structure`. Uses the k_dist = m - k_kofn + 1
conversion between kofn's :F convention (k = number of failures) and
dist.structure's :G convention (k = number of functioning components).

### Inclusion-exclusion: `R/ie_expand.R`
For exponential parallel systems, `prod(1 - exp(-lam_i * t))` expands
into a signed sum of exponentials. This makes all integrals closed-form.
O(2^m) terms, practical for m <= ~15.

### Observation schemes: `R/observe.R`
Ordered from least to most informative:
- Scheme 0: `observe_right_censor()`, system lifetime only (black box)
- Scheme 1: `observe_periodic(delta)`, periodic inspection
- Scheme 2: `observe_exact()`, complete monitoring (trivial)

### Estimation methods
- **Exponential** (`R/exp_kofn.R`): IE expansion for parallel, general
  critical-state enumeration for arbitrary k. Direct MLE via L-BFGS-B
  with Nelder-Mead fallback.
- **Weibull** (`R/wei_kofn.R`): EM algorithm (E-step via incomplete
  gamma, M-step via profile optimization) or direct MLE. Truncated
  moments in `R/truncated_moments.R`.
- **Scheme 1** (`R/scheme1.R`): Periodic inspection likelihood for
  both exponential and Weibull. Resolves shape-scale ambiguity.

### Data conventions
- `t`: system lifetime column
- `omega`: observation type ("exact", "right", "left", "interval")
- `x1, x2, ..., xm`: boolean candidate set columns
- `t_upper`: interval upper bound (for interval censoring)
- `comp_lower_j, comp_upper_j`: component inspection intervals (Scheme 1)

### Parameter conventions
- Exponential: `par = c(rate_1, ..., rate_m)`, m rates
- Weibull: `par = c(shape_1, scale_1, ..., shape_m, scale_m)`, 2m params interleaved
- All parameters must be positive
- Ascending parameter ordering for identifiability under permutation symmetry

### Component family dispatch
The `component` argument to `kofn()` is a `dfr_dist` prototype from
flexhaz. Family-specific S3 subclasses (`"dfr_exponential"`,
`"dfr_weibull"`) drive dispatch in `kofn_subclass()`,
`parse_params()`, `n_par_kofn()`, and the internal topology helpers.
Adding a new supported family requires:
1. One entry in `kofn_subclass()` mapping the dfr subclass to a kofn subclass
2. S3 methods for loglik, score, hess_loglik, fit, rdata on the new kofn subclass

## Key Mathematical Context

### k-out-of-n as censoring
For a k-out-of-n system with lifetime T = T_{(k)} (k-th order statistic):
- 1 component observed exactly (the critical component)
- k-1 components left-censored at T (already failed)
- n-k components right-censored at T (still functioning)

### Exponential parallel (k=m)
- f_sys(t) = sum_j lambda_j * exp(-lambda_j * t) * prod_{i!=j} (1 - exp(-lambda_i * t))
- IE expansion makes prod(1-exp) a finite sum of exponentials
- All integrals closed-form -> analytical likelihood

### Weibull parallel (k=m)
- f_sys(t) = sum_j f_j(t) * prod_{i!=j} F_i(t)
- No closed-form IE expansion
- EM: J = argmax_j T_j is latent
- E-step: truncated Weibull moments via incomplete gamma
- M-step: profile over shape, closed-form scale

## Testing
- Test exponential IE expansion against dist.structure system density
- Cross-validate with maskedcauses for the series case (k=1)
- Monte Carlo parameter recovery tests (stress tests, `skip_on_cran()`)
- Convergence tests for EM algorithm
