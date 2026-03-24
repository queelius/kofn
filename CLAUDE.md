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

- `maskedcauses`: series systems (k=1), masked cause of failure
- `maskedhaz`: general DFR framework for series systems
- `kofn`: parallel (k=m) and general k-out-of-n systems
- `likelihood.model`: generics (loglik, score, hess_loglik, fit, rdata)
- `generics`: fit generic
- `serieshaz` / `flexhaz`: distribution infrastructure

## Commands

```bash
Rscript -e "devtools::load_all()"          # Interactive development
Rscript -e "devtools::document()"          # Generate roxygen2 docs
Rscript -e "devtools::test()"              # Run tests
Rscript -e "covr::package_coverage()"      # Coverage report
Rscript -e "devtools::check()"             # R CMD check
```

## Architecture

### System representation: `R/coherent_system.R`
Coherent systems via minimal path sets. Precomputes minimal cut sets via
Berge transversal algorithm. Standard constructors: `parallel_system()`,
`series_system()`, `kofn_system()`, `bridge_system()`.

Key function: `system_censoring(system, times)`. Given component
lifetimes, returns system lifetime plus per-component censoring status
(exact/left/right).

### Inclusion-exclusion: `R/ie_expand.R`
For exponential parallel systems, `prod(1 - exp(-lam_i * t))` expands
into a signed sum of exponentials. This makes all integrals closed-form.
O(2^m) terms, practical for m <= ~15.

### Model API: `R/kofn.R`
Single constructor `kofn(k, m, family, method)` returns an S3 object.
Classes: `exp_kofn` (exponential) or `wei_kofn` (Weibull), inheriting
from `kofn` and `likelihood_model`.

Generic methods (from likelihood.model) return closures:
- `loglik(model)` -> `function(df, par)`
- `score(model)` -> `function(df, par)`
- `fit(model)` -> `function(df, par0, ...)`
- `rdata(model)` -> `function(theta, n, ...)`

### Observation schemes: `R/observe.R`
Ordered from least to most informative:
- Scheme 0: `observe_right_censor()`, system lifetime only (black box)
- Scheme 1: `observe_periodic(delta)`, periodic inspection
- Scheme 2: `observe_exact()`, complete monitoring (trivial)

### Estimation methods
- **Exponential** (`R/exp_kofn.R`): IE expansion for parallel, general
  critical-state enumeration for arbitrary k. Direct MLE via L-BFGS-B
  with Nelder-Mead fallback. Also has exponential EM for Scheme 0.
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
- Test exponential IE expansion against general system density engine
- Cross-validate with maskedcauses for the series case (k=1)
- Monte Carlo parameter recovery tests
- Convergence tests for EM algorithm
