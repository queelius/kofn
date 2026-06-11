# Package index

## Model constructor

Create a likelihood model for k-out-of-n systems.

- [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md) : Create
  a k-out-of-n system estimation model
- [`is_kofn()`](https://queelius.github.io/kofn/reference/is_kofn.md) :
  Test if an object is a kofn model
- [`ncomponents(`*`<kofn>`*`)`](https://queelius.github.io/kofn/reference/ncomponents.kofn.md)
  : Number of components in a kofn model
- [`print(`*`<kofn>`*`)`](https://queelius.github.io/kofn/reference/print.kofn.md)
  : Print method for kofn models

## Likelihood generics

S3 methods for the likelihood.model protocol.

- [`loglik(`*`<exp_kofn>`*`)`](https://queelius.github.io/kofn/reference/loglik.exp_kofn.md)
  : Log-likelihood for exponential k-out-of-n model
- [`loglik(`*`<wei_kofn>`*`)`](https://queelius.github.io/kofn/reference/loglik.wei_kofn.md)
  : Log-likelihood for Weibull k-out-of-n system
- [`score(`*`<exp_kofn>`*`)`](https://queelius.github.io/kofn/reference/score.exp_kofn.md)
  : Score function for exponential k-out-of-n model
- [`score(`*`<wei_kofn>`*`)`](https://queelius.github.io/kofn/reference/score.wei_kofn.md)
  : Score function for Weibull k-out-of-n system
- [`hess_loglik(`*`<exp_kofn>`*`)`](https://queelius.github.io/kofn/reference/hess_loglik.exp_kofn.md)
  : Hessian of the log-likelihood for exponential k-out-of-n model
- [`hess_loglik(`*`<wei_kofn>`*`)`](https://queelius.github.io/kofn/reference/hess_loglik.wei_kofn.md)
  : Hessian of log-likelihood for Weibull k-out-of-n system
- [`fit(`*`<exp_kofn>`*`)`](https://queelius.github.io/kofn/reference/fit.exp_kofn.md)
  : MLE fitting for exponential k-out-of-n model
- [`fit(`*`<wei_kofn>`*`)`](https://queelius.github.io/kofn/reference/fit.wei_kofn.md)
  : Fit Weibull k-out-of-n system model
- [`rdata(`*`<exp_kofn>`*`)`](https://queelius.github.io/kofn/reference/rdata.exp_kofn.md)
  : Random data generation for exponential k-out-of-n model
- [`rdata(`*`<wei_kofn>`*`)`](https://queelius.github.io/kofn/reference/rdata.wei_kofn.md)
  : Generate random data from a Weibull k-out-of-n system
- [`assumptions(`*`<exp_kofn>`*`)`](https://queelius.github.io/kofn/reference/assumptions.exp_kofn.md)
  : Assumptions for exponential k-out-of-n model
- [`assumptions(`*`<wei_kofn>`*`)`](https://queelius.github.io/kofn/reference/assumptions.wei_kofn.md)
  : Assumptions for Weibull k-out-of-n system model

## Observation functors

Composable observation schemes for censored and inspected data.

- [`observe_exact()`](https://queelius.github.io/kofn/reference/observe_exact.md)
  : Exact observation scheme (no censoring)
- [`observe_right_censor()`](https://queelius.github.io/kofn/reference/observe_right_censor.md)
  : Right-censoring observation scheme
- [`observe_left_censor()`](https://queelius.github.io/kofn/reference/observe_left_censor.md)
  : Left-censoring observation scheme
- [`observe_interval_censor()`](https://queelius.github.io/kofn/reference/observe_interval_censor.md)
  : Interval-censoring observation scheme
- [`observe_periodic()`](https://queelius.github.io/kofn/reference/observe_periodic.md)
  : Periodic inspection observation scheme
- [`observe_mixture()`](https://queelius.github.io/kofn/reference/observe_mixture.md)
  : Mixture of observation schemes

## Periodic inspection (Scheme 1)

- [`loglik_scheme1()`](https://queelius.github.io/kofn/reference/loglik_scheme1.md)
  : Log-likelihood for Scheme 1 (periodic inspection) parallel system
- [`rdata_scheme1()`](https://queelius.github.io/kofn/reference/rdata_scheme1.md)
  : Generate Scheme 1 (periodic inspection) data
- [`fit_scheme1()`](https://queelius.github.io/kofn/reference/fit_scheme1.md)
  : Fit k-out-of-n system MLE under Scheme 1 (periodic inspection)

## Masked cause-of-failure

- [`loglik_masked()`](https://queelius.github.io/kofn/reference/loglik_masked.md)
  : Masked k-out-of-n log-likelihood
- [`rdata_masked()`](https://queelius.github.io/kofn/reference/rdata_masked.md)
  : Generate masked k-out-of-n data

## Fisher information

- [`compare_fisher_info()`](https://queelius.github.io/kofn/reference/compare_fisher_info.md)
  : Compare Fisher information across observation schemes

## Parallel fast-path internals

Inclusion-exclusion expansion for parallel exponential systems.

- [`ie_expand()`](https://queelius.github.io/kofn/reference/ie_expand.md)
  : Inclusion-exclusion expansion of a product of CDFs
- [`w_j_exact()`](https://queelius.github.io/kofn/reference/w_j_exact.md)
  : Compute w_j(t) = f_j(t) \* prod_i != j F_i(t) for exponential
  components
- [`w_j_integral()`](https://queelius.github.io/kofn/reference/w_j_integral.md)
  : Closed-form integral of w_j(t) over an interval
- [`F_sys_exp()`](https://queelius.github.io/kofn/reference/F_sys_exp.md)
  : System CDF for exponential parallel systems
- [`S_sys_exp()`](https://queelius.github.io/kofn/reference/S_sys_exp.md)
  : System survival function for exponential parallel systems

## Shared helpers

- [`solve_mle()`](https://queelius.github.io/kofn/reference/solve_mle.md)
  : Default MLE solver for positive parameters
- [`trunc_log_moment_vec()`](https://queelius.github.io/kofn/reference/trunc_log_moment_vec.md)
  : Vectorized truncated log-moment of the Weibull distribution
- [`trunc_pow_moment()`](https://queelius.github.io/kofn/reference/trunc_pow_moment.md)
  : Scalar truncated power moment of the Weibull distribution
- [`trunc_pow_moment_vec()`](https://queelius.github.io/kofn/reference/trunc_pow_moment_vec.md)
  : Vectorized truncated power moment of the Weibull distribution
