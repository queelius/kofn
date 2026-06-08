## R CMD check results

0 errors | 0 warnings | 1 note

* One NOTE: "unable to verify current time" (system clock check, benign).

## First CRAN submission (v0.4.0)

kofn provides maximum likelihood estimation for k-out-of-n system data
in the Fisherian tradition. A k-out-of-n system fails when k of m
components have failed; the system lifetime is the k-th order statistic
of the component lifetimes, which creates structured censoring on
individual components (one observed exactly, k-1 left-censored, m-k
right-censored). The package supports exponential and Weibull components
under three observation schemes (system lifetime only, periodic
inspection, complete monitoring), with an EM algorithm for the Weibull
parallel case where direct MLE is unreliable.

The package conforms to the 'likelihood.model' generics and returns
fitted objects compatible with 'algebraic.mle'. The data-generating
process and topology infrastructure are delegated to 'dist.structure';
kofn focuses exclusively on inference for the k-out-of-n family.

All dependencies are already on CRAN:

* algebraic.dist (>= 1.0.0)
* algebraic.mle (>= 2.0.2)
* compositional.mle (>= 2.0.0)
* dist.structure (>= 0.5.0)
* flexhaz (>= 0.5.2)
* generics
* likelihood.model (>= 1.0.1)
* numDeriv

## Vignette runtime

Vignettes are configured to build under CRAN's per-vignette time
budget. Heavier Monte Carlo simulations are gated by a `run_long`
flag (default `FALSE`) so that the unconditional chunks complete in
roughly 30 to 100 seconds each; the full 7-vignette build is about
5 minutes locally.

## Test runtime

Local `devtools::test()` runs the full suite including stress tests
in roughly 8 minutes. The stress tests (`test-stress.R`,
`test-stress2.R`) use `skip_on_cran()` so that CRAN's incoming pretest
runs only the fast tests, totaling roughly 50 seconds.

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder (R-devel and R-release)
