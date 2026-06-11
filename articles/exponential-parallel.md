# Exponential Parallel Systems: Closed-Form MLE via Inclusion-Exclusion

## The Exponential Parallel System

A **parallel system** of $m$ independent components functions as long as
at least one component survives. The system lifetime is
$$T = \max\left( T_{1},\ldots,T_{m} \right),\qquad T_{j} \sim {Exp}\left( \lambda_{j} \right).$$

The system density is the derivative of the system CDF
$F_{sys}(t) = \prod_{j = 1}^{m}F_{j}(t) = \prod_{j = 1}^{m}\left( 1 - e^{- \lambda_{j}t} \right)$:
$$f_{sys}(t) = \sum\limits_{j = 1}^{m}f_{j}(t)\prod\limits_{i \neq j}F_{i}(t) = \sum\limits_{j = 1}^{m}\lambda_{j}e^{- \lambda_{j}t}\prod\limits_{i \neq j}(1 - e^{- \lambda_{i}t}).$$

Each summand $w_{j}(t) = f_{j}(t)\prod_{i \neq j}F_{i}(t)$ is the
contribution from the event “component $j$ is the last to fail.”

``` r
# Create an exponential parallel system model with 3 components
model <- kofn(k = 3, m = 3, component = dfr_exponential())
print(model)
#> k-out-of-n System Likelihood Model
#> -----------------------------------
#>   System type: parallel (k=3, m=3)
#>   Component distribution: exponential 
#>   Estimation method: mle 
#>   Column conventions:
#>     lifetime: t 
#>     omega: omega 
#>     interval upper: t_upper
```

The [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md)
constructor with `k = m` creates a parallel system. The model object
carries the system structure, column name conventions, and dispatches to
methods for
[`loglik()`](https://queelius.github.io/likelihood.model/reference/loglik.html),
[`score()`](https://queelius.github.io/likelihood.model/reference/score.html),
[`hess_loglik()`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html),
[`fit()`](https://generics.r-lib.org/reference/fit.html), and
[`rdata()`](https://queelius.github.io/likelihood.model/reference/rdata.html).

## Inclusion-Exclusion Expansion

The product $\prod_{i}\left( 1 - e^{- \lambda_{i}t} \right)$ cannot be
differentiated or integrated term-by-term as written. The
**inclusion-exclusion expansion** rewrites it as a finite signed sum of
exponentials:
$$\prod\limits_{i = 1}^{m}(1 - e^{- \lambda_{i}t}) = \sum\limits_{S \subseteq \{ 1,\ldots,m\}}( - 1)^{|S|}\, e^{- {(\sum\limits_{i \in S}\lambda_{i})}\, t}.$$

This has $2^{m}$ terms — practical for moderate $m$ (up to \$\$15
components). Because every term is a simple exponential $e^{- rt}$, all
integrals needed for censored-data likelihoods evaluate in closed form.

``` r
lam <- c(0.5, 0.3, 0.2)
ie <- ie_expand(lam)

# Display the expansion terms
data.frame(
  subset_size = sapply(seq_along(ie$sign), function(k) {
    sum(as.integer(intToBits(k - 1L))[1:3])
  }),
  sign = ie$sign,
  rate_sum = ie$rate_sum
)
#>   subset_size sign rate_sum
#> 1           0    1      0.0
#> 2           1   -1      0.5
#> 3           1   -1      0.3
#> 4           2    1      0.8
#> 5           1   -1      0.2
#> 6           2    1      0.7
#> 7           2    1      0.5
#> 8           3   -1      1.0
```

The first row ($S = \varnothing$) has sign $+ 1$ and rate sum $0$,
contributing the constant term $1$. The three singleton subsets
contribute $- e^{- \lambda_{j}t}$, the three pairs contribute
$+ e^{- {(\lambda_{i} + \lambda_{j})}t}$, and the full set contributes
$- e^{- {(\lambda_{1} + \lambda_{2} + \lambda_{3})}t}$.

We can verify the expansion numerically:

``` r
t_val <- 2.0

# Direct product
direct <- prod(1 - exp(-lam * t_val))

# Via IE expansion
ie_val <- sum(ie$sign * exp(-ie$rate_sum * t_val))

c(direct = direct, ie_expansion = ie_val, difference = direct - ie_val)
#>       direct ie_expansion   difference 
#>      0.09403      0.09403      0.00000
```

The key functions built on this expansion are:

- `w_j_exact(t, par, j)` — evaluates $w_{j}(t)$ at a single time point
- `w_j_integral(a, b, par, j)` — closed-form $\int_{a}^{b}w_{j}(t)\, dt$
- `F_sys_exp(t, par)` / `S_sys_exp(t, par)` — system CDF and survival

``` r
# The three w_j contributions at t = 2
sapply(1:3, function(j) w_j_exact(t_val, lam, j))
#> [1] 0.02736 0.03431 0.03824

# Their sum equals the system density
sum(sapply(1:3, function(j) w_j_exact(t_val, lam, j)))
#> [1] 0.09991

# System CDF and survival at t = 2
c(F_sys = F_sys_exp(t_val, lam), S_sys = S_sys_exp(t_val, lam))
#>   F_sys   S_sys 
#> 0.09403 0.90597
```

## Likelihood Under Different Observation Types

The IE expansion gives closed-form log-likelihood contributions for all
four observation types. Let
${\mathbf{λ}} = \left( \lambda_{1},\ldots,\lambda_{m} \right)$ and
$C_{i} \subseteq \{ 1,\ldots,m\}$ be the candidate set for observation
$i$.

**Exact** observation at $t$:
$$\ell_{i} = \log\sum\limits_{j \in C_{i}}w_{j}(t;\,{\mathbf{λ}}).$$

**Right-censored** at $t$ (system still functioning):
$$\ell_{i} = \log S_{sys}(t;\,{\mathbf{λ}}).$$

**Left-censored** at $t$ (system already failed):
$$\ell_{i} = \log\frac{\sum\limits_{j \in C_{i}}\int_{0}^{t}w_{j}(s)\, ds}{F_{sys}(t;\,{\mathbf{λ}})}.$$

**Interval-censored** on $(a,b\rbrack$:
$$\ell_{i} = \log\frac{\sum\limits_{j \in C_{i}}\int_{a}^{b}w_{j}(s)\, ds}{F_{sys}(b) - F_{sys}(a)}.$$

All integrals $\int_{a}^{b}w_{j}(s)\, ds$ are computed analytically by
[`w_j_integral()`](https://queelius.github.io/kofn/reference/w_j_integral.md).
The `loglik(model)` closure handles all four types automatically:

``` r
set.seed(42)

# Generate data under Scheme 0 (system lifetime only, no censoring)
gen <- rdata(model)
df <- gen(theta = lam, n = 50)
head(df, 5)
#>        t omega
#> 1  6.355 exact
#> 2  5.498 exact
#> 3  3.764 exact
#> 4  9.446 exact
#> 5 17.799 exact

# Evaluate the log-likelihood at the true parameters
ll_fn <- loglik(model)
ll_fn(df, lam)
#> [1] -144.9

# Compare with a wrong parameter vector
ll_fn(df, c(1, 1, 1))
#> [1] -320.3
```

## Maximum Likelihood Estimation

### Data generation

We generate $n = 100$ system lifetimes from a 3-component parallel
system with rates ${\mathbf{λ}} = (0.5,0.3,0.2)$ under Scheme 0 (system
lifetime only — no component-level information).

``` r
set.seed(123)
theta_true <- c(0.5, 0.3, 0.2)
n <- 50

gen <- rdata(model)
df <- gen(theta = theta_true, n = n)

cat("Observations:", nrow(df), "\n")
#> Observations: 50
cat("Observation types:\n")
#> Observation types:
table(df$omega)
#> 
#> exact 
#>    50
```

Under the default observation scheme (exact), we observe only when the
system fails, not which component failed last. This is the Scheme 0
(black-box) setting.

### Fitting

The `fit(model)` closure performs multi-start L-BFGS-B optimization with
Nelder-Mead fallback. Standard errors come from the observed Fisher
information (inverse of the negative Hessian at the MLE).

``` r
fit_fn <- fit(model)
result <- fit_fn(df, n_starts = 1L)

# MLE estimates
coef(result)
#> [1] 0.2806 0.2804 0.3752

# Standard errors
sqrt(diag(vcov(result)))
#> [1] 0.2911 0.2897 0.3313

# 95% confidence intervals
confint(result)
#>           2.5%  97.5%
#> param1 -0.2899 0.8511
#> param2 -0.2875 0.8482
#> param3 -0.2741 1.0244

# Log-likelihood, AIC, BIC
c(loglik = logLik(result), AIC = AIC(result), BIC = BIC(result))
#> loglik    AIC    BIC 
#> -128.4  262.9  268.6
```

``` r
summary(result)
#> Maximum likelihood estimator of type solver_result is normally distributed.
#> The estimates of the parameters are given by:
#> [1] 0.2806 0.2804 0.3752
#> The standard error is  0.2911 0.2897 0.3313 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>           2.5%  97.5%
#> param1 -0.2899 0.8511
#> param2 -0.2875 0.8482
#> param3 -0.2741 1.0244
#> The MSE of the individual components in a multivariate estimator is:
#>          [,1]     [,2]     [,3]
#> [1,]  0.08473 -0.07158 -0.02457
#> [2,] -0.07158  0.08395 -0.02359
#> [3,] -0.02457 -0.02359  0.10974
#> The log-likelihood is  -128.4 .
#> The AIC is  262.9 .
```

Every method used above, `coef`, `vcov`, `confint`, `logLik`, `AIC`,
`BIC`, `summary`, comes from `algebraic.mle` via its `mle` result class,
not from `kofn` itself. `kofn` supplies the likelihood; the inference
operations are inherited. See
[`vignette("ecosystem")`](https://queelius.github.io/kofn/articles/ecosystem.md)
for the full walk through the rlang MLE stack.

### Comparison with true values

``` r
est <- coef(result)
se <- sqrt(diag(vcov(result)))
ci <- confint(result)

comparison <- data.frame(
  component = 1:3,
  true = theta_true,
  estimate = est,
  se = se,
  ci_lower = ci[, 1],
  ci_upper = ci[, 2],
  covered = theta_true >= ci[, 1] & theta_true <= ci[, 2]
)
comparison
#>        component true estimate     se ci_lower ci_upper covered
#> param1         1  0.5   0.2806 0.2911  -0.2899   0.8511    TRUE
#> param2         2  0.3   0.2804 0.2897  -0.2875   0.8482    TRUE
#> param3         3  0.2   0.3752 0.3313  -0.2741   1.0244    TRUE
```

## Permutation Symmetry

Under Scheme 0, all components are candidates for every observation. The
log-likelihood is then **invariant under permutation** of the parameter
vector: swapping $\left. \lambda_{1}\leftrightarrow\lambda_{2} \right.$
yields the same log-likelihood value. This means there are $m!$
equivalent global maxima.

``` r
ll_fn <- loglik(model)

# Original parameter order
ll_original <- ll_fn(df, theta_true)

# Permuted parameters: swap components 1 and 3
theta_permuted <- theta_true[c(3, 1, 2)]
ll_permuted <- ll_fn(df, theta_permuted)

# Reversed order
theta_reversed <- rev(theta_true)
ll_reversed <- ll_fn(df, theta_reversed)

data.frame(
  ordering = c("original (0.5, 0.3, 0.2)",
               "permuted (0.2, 0.5, 0.3)",
               "reversed (0.2, 0.3, 0.5)"),
  loglik = c(ll_original, ll_permuted, ll_reversed),
  difference = c(0, ll_permuted - ll_original, ll_reversed - ll_original)
)
#>                   ordering loglik difference
#> 1 original (0.5, 0.3, 0.2) -129.1          0
#> 2 permuted (0.2, 0.5, 0.3) -129.1          0
#> 3 reversed (0.2, 0.3, 0.5) -129.1          0
```

The log-likelihood values are identical (up to machine precision). This
symmetry has two important consequences:

1.  **Initialization matters.** The optimizer may converge to any of the
    $m!$ equivalent modes. The
    [`default_init_exp()`](https://queelius.github.io/kofn/reference/default_init_exp.md)
    function uses a **spread initialization** — rates spaced as
    `lam0 * seq(0.5, 1.5, length.out = m)` — to break the symmetry and
    consistently find the ascending-order mode.

2.  **Sorted matching for Monte Carlo.** When evaluating estimator
    performance across replications, we sort both the estimates and the
    true parameters in ascending order before computing bias, RMSE, and
    coverage. This resolves the label-switching problem.

## Monte Carlo Validation

We assess the MLE’s finite-sample performance via Monte Carlo
simulation. For each replicate, we generate data from the true model,
fit the MLE, and record the estimates, standard errors, and confidence
interval coverage.

``` r
# Small demo: R = 5 replications for illustration (keeps build < 30 sec).
# Set run_long = TRUE above for the full R = 100 simulation.
set.seed(2026)

R <- 3
n_mc <- 50
alpha <- 0.05

theta_true_mc <- c(0.5, 0.3, 0.2)
theta_sorted <- sort(theta_true_mc)
m_mc <- length(theta_true_mc)

gen_mc <- rdata(model)
fit_mc <- fit(model)

estimates <- matrix(NA, nrow = R, ncol = m_mc)
se_hat <- matrix(NA, nrow = R, ncol = m_mc)
ci_lower <- matrix(NA, nrow = R, ncol = m_mc)
ci_upper <- matrix(NA, nrow = R, ncol = m_mc)
converged <- logical(R)

for (r in seq_len(R)) {
  df_r <- gen_mc(theta = theta_true_mc, n = n_mc)
  res_r <- tryCatch(fit_mc(df_r, n_starts = 1L), error = function(e) NULL)

  if (!is.null(res_r) && !any(is.na(coef(res_r)))) {
    ord <- order(coef(res_r))
    estimates[r, ] <- coef(res_r)[ord]
    se_r <- sqrt(diag(vcov(res_r)))
    se_hat[r, ] <- se_r[ord]
    ci_r <- confint(res_r, level = 1 - alpha)
    ci_lower[r, ] <- ci_r[ord, 1]
    ci_upper[r, ] <- ci_r[ord, 2]
    converged[r] <- TRUE
  }
}
```

### Results

After sorting each replicate’s estimates in ascending order (to match
the sorted true values
$\lambda_{(1)} = 0.2 < \lambda_{(2)} = 0.3 < \lambda_{(3)} = 0.5$), we
compute bias, RMSE, and coverage:

``` r
ok <- converged & complete.cases(estimates)
cat("Converged replications:", sum(ok), "of", R, "\n\n")
#> Converged replications: 3 of 3

est_ok <- estimates[ok, , drop = FALSE]
se_ok <- se_hat[ok, , drop = FALSE]
ci_lo <- ci_lower[ok, , drop = FALSE]
ci_hi <- ci_upper[ok, , drop = FALSE]

bias <- colMeans(est_ok) - theta_sorted
rmse <- sqrt(colMeans((est_ok - matrix(theta_sorted, nrow = sum(ok),
                                        ncol = m_mc, byrow = TRUE))^2))
coverage <- colMeans(
  ci_lo <= matrix(theta_sorted, nrow = sum(ok), ncol = m_mc, byrow = TRUE) &
  ci_hi >= matrix(theta_sorted, nrow = sum(ok), ncol = m_mc, byrow = TRUE)
)
mean_se <- colMeans(se_ok)

mc_table <- data.frame(
  component = paste0("lambda_(", 1:m_mc, ")"),
  true = theta_sorted,
  mean_est = colMeans(est_ok),
  bias = bias,
  rmse = rmse,
  mean_se = mean_se,
  coverage_95 = coverage
)
mc_table
#>    component true mean_est     bias    rmse mean_se coverage_95
#> 1 lambda_(1)  0.2   0.2151  0.01513 0.02827 0.06274      1.0000
#> 2 lambda_(2)  0.3   0.2888 -0.01124 0.08127 0.11225      1.0000
#> 3 lambda_(3)  0.5   0.6770  0.17699 0.54061 0.60901      0.6667
```

Key observations:

- **Bias** should be near zero for a correctly specified MLE.
- **RMSE** reflects finite-sample variability; the smallest rate
  ($\lambda_{(1)} = 0.2$) is typically estimated with larger relative
  error because it contributes least to the system density (its
  component is least likely to be the last to fail).
- **Coverage** near 0.95 validates the Hessian-based confidence
  intervals. Note that the SE ratio artifact — Hessian-based SEs
  computed at the unsorted MLE, then reordered — can cause mild under-
  or over-coverage.

## Observation Functors

The `observe_*` family of functions controls how system lifetimes are
recorded. Each functor maps a true system lifetime to an observed record
(time, observation type, and optional interval bounds).

### Scheme 0: System lifetime with right-censoring

`observe_right_censor(tau)` is the default: systems failing before
$\tau$ are observed exactly, and those surviving past $\tau$ are
right-censored.

``` r
# Different censoring times
obs_light <- observe_right_censor(tau = 20)   # Light censoring
obs_heavy <- observe_right_censor(tau = 5)    # Heavy censoring
obs_none  <- observe_right_censor()           # No censoring (tau = Inf)

# Example: system fails at t = 8
obs_light(8)   # exact
#> $t
#> [1] 8
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
obs_heavy(8)   # right-censored at tau = 5
#> $t
#> [1] 5
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA
obs_none(8)    # exact
#> $t
#> [1] 8
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
```

### Effect of censoring on estimation

Heavier right-censoring discards information about long-lived systems,
increasing estimator variance — especially for the smallest rate (the
most reliable component, which is most likely to be the last survivor).

``` r
set.seed(99)

tau_values <- c(5, Inf)

cat(sprintf("%-8s  %6s  %10s  %10s  %10s\n",
            "tau", "n_right", "est_lam1", "est_lam2", "est_lam3"))
#> tau       n_right    est_lam1    est_lam2    est_lam3
cat(strrep("-", 56), "\n")
#> --------------------------------------------------------

for (tau in tau_values) {
  obs_fn <- if (is.finite(tau)) observe_right_censor(tau) else NULL
  df_tau <- gen(theta = theta_true, n = 50,
                observe = obs_fn)
  n_right <- sum(df_tau$omega == "right")

  res_tau <- tryCatch(
    fit_fn(df_tau, n_starts = 1L),
    error = function(e) NULL
  )

  if (!is.null(res_tau)) {
    est <- sort(coef(res_tau))
    cat(sprintf("%-8s  %6d  %10.4f  %10.4f  %10.4f\n",
                if (is.finite(tau)) as.character(tau) else "Inf",
                n_right, est[1], est[2], est[3]))
  }
}
#> 5             28      0.1229      0.6073      3.3339
#> Inf            0      0.2554      0.2554      0.2554
```

### Scheme 1: Periodic inspection

Under `observe_periodic(delta, tau)`, the system is inspected at times
$\delta,2\delta,3\delta,\ldots$ The failure time is known only to lie
within the inspection interval containing it (interval-censored).

``` r
obs_periodic <- observe_periodic(delta = 2, tau = 30)

# System fails at t = 7.3 -> interval-censored to (6, 8]
obs_periodic(7.3)
#> $t
#> [1] 6
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 8

# System survives past tau = 30 -> right-censored
obs_periodic(35)
#> $t
#> [1] 30
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA
```

``` r
set.seed(42)
df_s1 <- gen(theta = theta_true, n = 50,
             observe = observe_periodic(delta = 2, tau = 50))
table(df_s1$omega)
#> 
#> interval 
#>       50

# Fit under interval censoring
res_s1 <- fit_fn(df_s1, n_starts = 1L)
sort(coef(res_s1))
#> [1] 0.1527 0.4575 0.5817
```

Interval censoring retains more information than pure right-censoring at
the same study duration, because every failure is localized to a finite
window rather than lost entirely.
