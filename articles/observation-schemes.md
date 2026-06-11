# Observation Scheme Composability

## The Observe Functor API

Every `observe_*` function returns a closure
`function(t_true) -> list(t, omega, t_upper)` that maps a true system
lifetime to an observed record. The `omega` field classifies the
observation: `"exact"`, `"right"`, `"left"`, or `"interval"`.

``` r
# Exact: the trivial case -- no information loss
observe_exact()(3.7)
#> $t
#> [1] 3.7
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA

# Right-censoring: systems surviving past tau are censored
obs_rc <- observe_right_censor(tau = 5)
obs_rc(3.7)   # fails before tau -> exact
#> $t
#> [1] 3.7
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
obs_rc(8.2)   # survives past tau -> right-censored at 5
#> $t
#> [1] 5
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA

# Left-censoring: systems failing before tau are censored
obs_lc <- observe_left_censor(tau = 2)
obs_lc(1.5)   # fails before tau -> left-censored at 2
#> $t
#> [1] 2
#> 
#> $omega
#> [1] "left"
#> 
#> $t_upper
#> [1] NA
obs_lc(3.7)   # fails after tau  -> exact
#> $t
#> [1] 3.7
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA

# Interval-censoring: failures in [a, b) are binned
obs_ic <- observe_interval_censor(a = 2, b = 6)
obs_ic(4.0)   # inside window  -> interval [2, 6)
#> $t
#> [1] 2
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 6
obs_ic(1.0)   # outside window -> exact
#> $t
#> [1] 1
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA

# Periodic inspection: regular grid with right-censoring at tau
obs_per <- observe_periodic(delta = 3, tau = 15)
obs_per(7.3)  # falls in [6, 9) -> interval
#> $t
#> [1] 6
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 9
obs_per(20)   # past tau        -> right-censored at 15
#> $t
#> [1] 15
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA
```

## Composing Schemes with rdata

The `rdata(model)` closure accepts any observe functor via its `observe`
argument. This decouples the data-generating process from the
observation protocol.

``` r
model <- kofn(k = 2, m = 2, component = dfr_exponential())
theta <- c(1.0, 0.5)
gen <- rdata(model)

# Exact observation (default)
df_exact <- gen(theta, n = 6)
head(df_exact)
#>        t omega
#> 1 0.6280 exact
#> 2 0.8203 exact
#> 3 2.3832 exact
#> 4 1.4297 exact
#> 5 2.6894 exact
#> 6 4.8174 exact

# Right-censoring at tau = 2
df_right <- gen(theta, n = 6, observe = observe_right_censor(tau = 2))
head(df_right)
#>        t omega
#> 1 2.0000 right
#> 2 0.7407 exact
#> 3 2.0000 right
#> 4 1.3317 exact
#> 5 2.0000 right
#> 6 0.6221 exact

# Periodic inspection every delta = 1 time unit
df_per <- gen(theta, n = 6, observe = observe_periodic(delta = 1, tau = 10))
head(df_per)
#>   t    omega t_upper
#> 1 1 interval       2
#> 2 6 interval       7
#> 3 1 interval       2
#> 4 0 interval       1
#> 5 3 interval       4
#> 6 1 interval       2
```

The same model and fitting infrastructure handles all observation types
transparently – the log-likelihood dispatches on the `omega` column.

## Effect on Estimation

We compare MLE quality across observation schemes for a 2-component
exponential parallel system with rates $\lambda = (1.0,0.5)$. We run R =
3 replicates of n = 60 observations under each scheme and report the
RMSE of the sorted rate estimates. This is a small demo for vignette
build speed; for tighter estimates set `R <- 100, n <- 200` and rerun.

``` r
set.seed(2026)
R <- 3; n <- 60
theta <- c(1.0, 0.5)
theta_sorted <- sort(theta)
model <- kofn(k = 2, m = 2, component = dfr_exponential())
gen <- rdata(model)
fit_fn <- fit(model)

schemes <- list(
  exact       = NULL,
  right_tau3  = observe_right_censor(tau = 3),
  right_tau1  = observe_right_censor(tau = 1),
  periodic_d1 = observe_periodic(delta = 1, tau = 20),
  left_tau1   = observe_left_censor(tau = 1)
)

results <- lapply(names(schemes), function(nm) {
  ests <- matrix(NA, nrow = R, ncol = 2)
  for (r in seq_len(R)) {
    df <- gen(theta, n, observe = schemes[[nm]])
    res <- tryCatch(fit_fn(df, n_starts = 1L), error = function(e) NULL)
    if (!is.null(res) && !any(is.na(coef(res))))
      ests[r, ] <- sort(coef(res))
  }
  ok <- complete.cases(ests)
  errs <- sweep(ests[ok, , drop = FALSE], 2, theta_sorted)
  data.frame(
    scheme    = nm,
    rmse_1    = round(sqrt(mean(errs[, 1]^2)), 3),
    rmse_2    = round(sqrt(mean(errs[, 2]^2)), 3),
    converged = sum(ok),
    stringsAsFactors = FALSE
  )
})
scheme_rmse <- do.call(rbind, results)
scheme_rmse
#>        scheme rmse_1 rmse_2 converged
#> 1       exact  0.130  0.372         3
#> 2  right_tau3  0.127  3.083         3
#> 3  right_tau1  0.220  2.369         3
#> 4 periodic_d1  0.040  4.670         3
#> 5   left_tau1  0.130  0.524         3
```

Key patterns:

- **Right-censoring** degrades as `tau` decreases – more observations
  are lost to censoring, reducing the effective sample size.
- **Periodic inspection** trades exact times for interval bounds. The
  interval-censored likelihood surface is flatter, which can cause
  occasional optimizer outliers that inflate the RMSE.
- **Left-censoring** at `tau = 1` has little effect here because most
  parallel system lifetimes exceed 1 (mean $\approx 2.3$).
  Left-censoring would hurt more if `tau` were closer to the median
  system lifetime.

Note that the parallel system’s permutation symmetry dominates these
differences – the optimizer’s ability to separate two rates depends more
on the likelihood geometry than on the observation scheme.

## Mixed Observation Environments

[`observe_mixture()`](https://queelius.github.io/kofn/reference/observe_mixture.md)
randomly selects a scheme for each observation. This models
heterogeneous monitoring – e.g., 70% of units are continuously
monitored, 30% are only inspected periodically.

``` r
obs_mix <- observe_mixture(
  observe_exact(),
  observe_periodic(delta = 2, tau = 20),
  weights = c(0.7, 0.3)
)

set.seed(99)
model <- kofn(k = 2, m = 2, component = dfr_exponential())
gen <- rdata(model)
df_mix <- gen(c(1.0, 0.5), n = 60, observe = obs_mix)
table(df_mix$omega)
#> 
#>    exact interval 
#>       42       18
```

The mixture preserves the functor interface – `fit(model)` handles the
resulting data identically, dispatching on the per-row `omega` values.

``` r
fit_fn <- fit(model)
res_mix <- fit_fn(df_mix, n_starts = 1L)
sort(coef(res_mix))
#> [1] 0.4797 7.1992
```
