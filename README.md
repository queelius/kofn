# kofn

**Component Lifetime Estimation from k-out-of-n System Data**

<!-- badges: start -->
<!-- badges: end -->

## The Problem

You observe a system — an engine, a power grid, a redundant server
cluster — and you record when it fails. But you don't care about the
*system*. You care about the *components* inside it. How long do the
individual bearings, transformers, or nodes last?

This is fundamentally a **censoring problem**. When a k-out-of-n system
fails (meaning k components have failed), the component lifetimes are
only partially observed:

- **1 component** is observed exactly (the one that triggered system
  failure)
- **k − 1 components** are right-censored (still functioning at system
  failure)
- **m − k components** are left-censored (already failed, but we don't
  know when)

The `kofn` package provides maximum likelihood estimation of component
lifetime parameters from this structured censoring, supporting
exponential and Weibull component distributions under multiple
observation schemes.

## Why This Matters

**Parallel systems** (k = n) are the hardest case: all components have
failed by the time the system fails, but we only observe the system
failure time — not the individual failure times. The last-to-fail
component is observed exactly; all others are deeply left-censored.

This creates an **information asymmetry**: fast components (those that
typically fail early) contribute almost no information to the likelihood,
because their CDF is near 1 at typical system failure times. The `kofn`
package quantifies this asymmetry and provides tools to resolve it
through periodic inspection.

## Installation

```r
# Install from GitHub
remotes::install_github("queelius/kofn")
```

## Quick Start

```r
library(kofn)

# Create a parallel system model (k=m) with 3 exponential components
model <- kofn(k = 3, m = 3, family = "exponential")

# Generate data: true rates = (0.5, 0.3, 0.2)
gen <- rdata(model)
set.seed(42)
df <- gen(theta = c(0.5, 0.3, 0.2), n = 300)

# Fit via maximum likelihood
fitter <- fit(model)
result <- fitter(df, n_starts = 5)

# Full base R stats interface
coef(result)       # parameter estimates
confint(result)    # confidence intervals
summary(result)    # coefficient table with SEs
AIC(result)        # model comparison
```

## Weibull Components with EM

When components follow a Weibull distribution (capturing increasing or
decreasing failure rates), the estimation problem is substantially
harder. The `kofn` package provides an EM algorithm whose E-step
computes truncated Weibull moments via incomplete gamma functions.

```r
# Weibull parallel system with EM estimation
model_wei <- kofn(k = 2, m = 2, family = "weibull", method = "em")
gen <- rdata(model_wei)
set.seed(42)
df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 300)
#                   ^shape ^scale ^shape ^scale

fitter <- fit(model_wei)
result <- fitter(df, n_starts = 3)
summary(result)
```

## Periodic Inspection (Scheme 1)

Under Scheme 0 (system-level observation only), fast components are
poorly estimated. Periodic inspection at interval $\delta$ localizes
each component's failure time to an inspection interval, resolving the
information asymmetry:

```r
# Generate data with periodic inspection every 0.5 time units
s1_gen <- rdata_scheme1(model_wei)
df_s1 <- s1_gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 300, delta = 0.5)

# Fit using the Scheme 1 likelihood
s1_fitter <- fit_scheme1(model_wei)
result_s1 <- s1_fitter(df_s1, n_starts = 3)
```

Shape RMSE drops by an order of magnitude — even coarse inspection
intervals provide dramatic improvement for the worst-estimated
components.

## Coherent System Engine

The package includes a general coherent system engine based on minimal
path and cut sets, supporting arbitrary system structures beyond
k-out-of-n:

```r
# Standard systems
sys_parallel <- parallel_system(5)
sys_series   <- series_system(5)
sys_2of5     <- kofn_system(2, 5)
sys_bridge   <- bridge_system()

# System properties
system_signature(sys_2of5)
phi(sys_bridge, c(TRUE, FALSE, TRUE, TRUE, FALSE))
system_censoring(sys_parallel, c(1.2, 3.4, 2.1, 0.8, 4.5))
```

## Key Features

| Feature | Description |
|---------|-------------|
| **Exponential MLE** | Closed-form likelihood via inclusion-exclusion expansion |
| **Weibull EM** | EM algorithm with incomplete gamma E-step and profile M-step |
| **Direct MLE** | L-BFGS-B with Nelder-Mead fallback, multi-start |
| **Observation schemes** | Scheme 0 (black box), Scheme 1 (periodic inspection), Scheme 2 (complete) |
| **Fisher information** | Compare information across observation schemes |
| **Coherent systems** | Arbitrary systems via minimal path/cut sets |
| **Base R integration** | `coef()`, `vcov()`, `logLik()`, `AIC()`, `confint()`, `summary()` |

## Ecosystem

`kofn` is part of a family of packages for likelihood-based reliability
inference:

- [**maskedcauses**](https://github.com/queelius/maskedcauses) — Series
  systems (k = n) with masked cause of failure
- [**likelihood.model**](https://github.com/queelius/likelihood.model) —
  Likelihood-based inference generics and `fisher_mle` objects
- [**algebraic.mle**](https://github.com/queelius/algebraic.mle) — MLE
  objects with algebraic operations

## Vignettes

- `vignette("censoring-framework")` — The theoretical foundation:
  k-out-of-n as structured censoring
- `vignette("exponential-parallel")` — Exponential parallel systems with
  inclusion-exclusion
- `vignette("kofn-systems")` — Coherent systems and the k-out-of-n
  spectrum
- `vignette("weibull-em")` — Weibull EM algorithm and why exponential
  theory doesn't generalize
- `vignette("periodic-inspection")` — Observation schemes and resolving
  information asymmetry
