# Inside kofn: Building on the rlang MLE Stack

## The design in one paragraph

`kofn` is a small package. Most of what a user sees when working with
it, the inference methods, the optimizer composition, the ability to
combine independent fits, the bridge to distribution algebra, is not
defined inside `kofn`. It comes from four sibling packages that sit
underneath `kofn` in a deliberately layered stack. This vignette is a
guided tour of that stack using a `kofn` model as the running example.
The point is to show what `kofn` had to specialize (the domain-specific
content) and what it gets for free by conforming to ecosystem-wide
generics.

The layering is a design choice, not an accident. Building reliability
packages on a uniform substrate means the next package up the chain,
`dagsys` for general coherent systems, or a future `maskedparallel`, can
reuse the same inference machinery without rebuilding it. Specialization
only enters when a new concept genuinely cannot be expressed in the
existing vocabulary.

## The stack

      kofn                 (this package: k-out-of-m censoring MLE)
        |
      likelihood.model     (likelihood calculus as a set of generics)
      compositional.mle    (MLE solvers as first-class composable objects)
        |
      algebraic.mle        (MLE objects with asymptotic inference)
        |
      algebraic.dist       (probability distribution algebra)

Each layer crystallizes a single concept:

| Layer               | Concept                                                                                                                                                                                                                                                    |
|---------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `algebraic.dist`    | A probability distribution is a first-class value. `normal(0, 1) + normal(0, 1)` is a `normal(0, 2)`. Distribution arithmetic simplifies when it can and falls back to Monte Carlo when it cannot.                                                         |
| `algebraic.mle`     | An MLE is not a number. It is an object with an asymptotic covariance. `confint`, `vcov`, `AIC`, `bias`, `mse`, `pred`, and `combine` are defined on it.                                                                                                   |
| `likelihood.model`  | A likelihood model is specified by generics (`loglik`, `score`, `hess_loglik`, `fim`, `rdata`, `assumptions`). A package that implements these generics gets Fisherian inference (support functions, relative likelihoods, likelihood intervals) for free. |
| `compositional.mle` | A solver is a function. Functions compose. `lbfgsb() %>>% nelder_mead()` is a sequential chain; `lbfgsb() %|% nelder_mead()` is a parallel race. `with_restarts(lbfgsb(), n = 5)` is a wrapper.                                                            |

The `compositional.mle` package has its own vignette, `mle-ecosystem`,
that walks through the layers from the other direction. Read it if you
want the bottom-up view. This vignette is top-down: start with a `kofn`
model and peel back the layers.

## A running example

We will use a 3-component series system. The system fails when the first
component fails (series, $k = 1$); component lifetimes are i.i.d.
exponential with unknown rates. Series is the cleanest case for
demonstrating the ecosystem generics; richer configurations (and the
identifiability issues that come with them) are the subject of other
vignettes.

``` r
library(kofn)
library(flexhaz)
set.seed(42)

model <- kofn(k = 1, m = 3, component = dfr_exponential())
gen   <- rdata(model)
df    <- gen(theta = c(0.3, 0.8, 1.5), n = 100)
head(df, 3)
#>            t omega
#> 1 0.54158264 exact
#> 2 0.08642963 exact
#> 3 0.16561603 exact
```

The object `model` is our entry point into the ecosystem. It is an S3
object inheriting from both `kofn` and `likelihood_model`.

``` r
class(model)
#> [1] "exp_kofn"         "kofn"             "likelihood_model"
```

That second class, `likelihood_model`, is the contract: `kofn` has
promised to implement the likelihood calculus as generics.

## Layer 1: the likelihood calculus (`likelihood.model`)

`likelihood.model` exports four small generics: `loglik`, `score`,
`hess_loglik`, and `fim`. Each one is a *closure factory*: you pass a
model and get back a function that takes data.

``` r
ll <- loglik(model)                  # ll: function(df, par) -> numeric
sc <- score(model)                   # sc: function(df, par) -> gradient
H  <- hess_loglik(model)             # H:  function(df, par) -> Hessian

par0 <- c(0.3, 0.8, 1.5)
ll(df, par0)
#> [1] -22.51077
sc(df, par0)
#> [1] -6.946891 -6.946891 -6.946891
```

The package implements `loglik.exp_kofn` and `score.exp_kofn` as S3
methods. The generic dispatch is what makes it possible for any other
package, one that only knows about `likelihood_model` objects, to
consume a `kofn` model without knowing anything about k-out-of-m
systems. This is the payoff of conforming to a published interface.

`likelihood.model` also provides `rdata` (a closure factory for
simulating data from the model) and `assumptions` (a description of what
the model assumes). `kofn` implements both.

``` r
assumptions(model)[1:3]              # first 3 assumptions
#> [1] "independent component lifetimes"             
#> [2] "exponential component lifetime distributions"
#> [3] "1-out-of-3 system structure"
```

## Layer 2: the optimizer (`compositional.mle`)

`fit(model)` returns a closure that optimizes the likelihood. Internally
it delegates to
[`kofn::solve_mle`](https://queelius.github.io/kofn/reference/solve_mle.md),
which wraps `compositional.mle`’s L-BFGS-B primary solver with a
log-scale Nelder-Mead fallback. That detail is intentionally opaque to
the caller.

``` r
fit_fn <- fit(model)
res <- fit_fn(df, n_starts = 1L)
round(coef(res), 3)
#> [1] 1.458 0.581 0.163
```

What if you want a different optimization strategy? You do not have to
write a new solver. You compose existing ones:

``` r
library(compositional.mle)
#> Loading required package: algebraic.mle

# Build the problem yourself to use compositional.mle directly
prob <- mle_problem(
  loglike    = function(par) ll(df, par),
  constraint = mle_constraint(support = function(par) all(par > 0))
)

# Sequential: grid search as a warm start, then L-BFGS-B for polish
strategy <- grid_search(lower = rep(0.05, 3),
                         upper = rep(2.0, 3), n = 3) %>>%
            lbfgsb(lower = rep(1e-10, 3))
res_chain <- strategy(prob, theta0 = rep(0.5, 3))
round(coef(res_chain), 3)
#>  Var1  Var2  Var3 
#> 0.084 0.084 2.034
```

The `%>>%` operator is **sequential composition**: run the first solver,
pipe its estimate into the second. There is also `%|%` (**parallel
race**), which runs two solvers in parallel and keeps whichever gets a
higher log-likelihood:

``` r
race <- lbfgsb(lower = rep(1e-10, 3)) %|% nelder_mead()
res_race <- race(prob, theta0 = rep(0.5, 3))
round(coef(res_race), 3)
#> [1] 0.734 0.734 0.734
```

Neither composition operator is defined in `kofn`. They apply to any MLE
problem specified through `mle_problem`. `kofn` simply built its
likelihood to be compatible.

## Layer 3: the MLE object (`algebraic.mle`)

Fit results carry inference operations. All of the following come from
`algebraic.mle` (or methods that `algebraic.mle` registers on its result
classes), not from `kofn`:

``` r
library(algebraic.mle)

round(coef(res), 3)                  # point estimates
#> [1] 1.458 0.581 0.163
round(se(res), 3)                    # standard errors
#> [1]  6.232 14.575 30.840
round(confint(res), 3)[1:3, ]        # 95% CIs, first three params
#>           2.5%  97.5%
#> param1 -10.757 13.673
#> param2 -27.985 29.146
#> param3 -60.282 60.608
logLik(res)                          # log-likelihood
#> 'log Lik.' -21.05276 (df=3)
AIC(res); BIC(res)                   # information criteria
#> [1] 48.10552
#> [1] 55.92103
nobs(res); nparams(res)              # bookkeeping
#> [1] 100
#> [1] 3
```

Beyond the classics, `algebraic.mle` supplies inference-theoretic
machinery:

``` r
round(bias(res), 4)                  # O(1/n) bias (zero for MLE at truth)
#> [1] 0 0 0
round(observed_fim(res), 1)[1:3, 1:3]  # observed Fisher information
#>      [,1] [,2] [,3]
#> [1,] 43.8 17.5  4.9
#> [2,] 17.5  7.0  2.0
#> [3,]  4.9  2.0  0.6
```

The payoff for bundling these under one class becomes clear when we
compose estimates. Suppose we run the same experiment twice (independent
replicates) and want to pool them:

``` r
df2  <- gen(theta = c(0.3, 0.8, 1.5), n = 100)
res1 <- fit_fn(df,  n_starts = 1L)
res2 <- fit_fn(df2, n_starts = 1L)

pooled <- combine(res1, res2)        # inverse-variance weighted
round(coef(pooled), 3)
#> [1]  2.207 -1.362  0.393
round(se(pooled), 3)
#> [1]  7.005  5.752 81.762
nobs(pooled)                         # 200 = 100 + 100
#> [1] 200
```

`combine` is not a `kofn` function. It is defined in `algebraic.mle` for
any MLE object whose asymptotic covariance is known. `kofn` gets it
because it uses the ecosystem’s currency.

## Layer 4: the distribution bridge (`algebraic.dist`)

MLE results are asymptotically normal.
[`algebraic.mle::as_dist`](https://queelius.github.io/algebraic.dist/reference/as_dist.html)
converts any MLE into the corresponding multivariate normal
distribution, at which point the full `algebraic.dist` algebra applies:

``` r
d <- as_dist(res)
class(d)
#> [1] "mvn"               "multivariate_dist" "continuous_dist"  
#> [4] "dist"
```

This matters when you want to reason about the sampling distribution of
derived quantities. For a redundant system, the mean time to the first
failure is a nonlinear function of the rates; its sampling distribution
can be obtained through the delta method, but
[`algebraic.mle::rmap`](https://queelius.github.io/algebraic.dist/reference/rmap.html)
does the bookkeeping automatically:

``` r
# g(lambda) = 1/sum(lambda): mean time to first failure (series system)
g <- function(lam) 1 / sum(lam)
res_mttff <- rmap(res, g)
round(coef(res_mttff), 3)            # point estimate
#> [1] 0.454
round(se(res_mttff), 4)              # delta-method SE
#> [1] 0.4461
```

The delta method is general. Nothing above is specific to `kofn`.

## What had to be specialized

Everything in the previous four sections is generic ecosystem code. So
what did `kofn` actually have to write? Three kinds of thing:

**1. The likelihood itself**, which encodes k-out-of-m censoring:

``` r
# R/exp_kofn.R (simplified)
loglik.exp_kofn <- function(model) {
  function(df, par) {
    # Inclusion-exclusion expansion for k = m (parallel).
    # Critical-state enumeration for general k.
    # Interval-censored contributions for Scheme 1 rows.
    ...
  }
}
```

This is where the domain work lives: the IE expansion for exponential
parallel, the critical-state enumeration for general k, the truncated
Weibull moments for the EM E-step, and the interval-censoring
contributions for Scheme 1. These could not be inherited because they
encode the physics of k-out-of-m censoring.

**2. Data-generating processes**. `rdata.exp_kofn` and `rdata.wei_kofn`
implement the forward simulator (generate component lifetimes, compute
system failure time, record censoring pattern). Trivial in principle,
but the package has to know the model.

**3. Custom observation schemes**. `loglik_scheme1` (periodic
inspection) and `loglik_masked` (partial autopsy) are two extra
observation mechanisms that require their own likelihood construction.
Each of them, once written, plugs back into the same generics above.

Everything else, the solver, the inference, the composition, the
distributional bridge, came from the ecosystem.

## Extending the stack

The same interface is what makes further specialization cheap. A future
`dagsys` package, for arbitrary coherent systems described as DAGs of
components, would do exactly what `kofn` does:

1.  Define a new S3 class (e.g., `dag_system`) inheriting from
    `likelihood_model`.
2.  Implement `loglik.dag_system`, `score.dag_system`, and
    `rdata.dag_system`.
3.  Arrange for `fit.dag_system` to call `solve_mle` (or a custom
    `compositional.mle` strategy).

Step 3 alone brings in `coef`, `vcov`, `confint`, `summary`, `logLik`,
`AIC`, `BIC`, `bias`, `mse`, `combine`, `rmap`, `as_dist`, and all of
distribution algebra. The specialization is only the domain-specific
likelihood, which is where the real work is anyway.

This is the bet the rlang stack is making: keep the generic layer small
and stable, and the specialization layer earns the right to exist by
introducing a genuinely new concept. `kofn` earned it by introducing
k-out-of-m censoring. `maskedcauses` earned it by introducing
series-system masked cause of failure and publishing the `series_md`
protocol (an S3 class specialising `likelihood_model`) that `maskedhaz`
then implements against a different hazard class. The shared
infrastructure is what lets both the protocol package and its
alternative implementations stay small.

## Further reading

- `vignette("mle-ecosystem", package = "compositional.mle")` walks the
  same stack bottom-up with simpler examples.
- [`vignette("getting-started", package = "likelihood.model")`](https://queelius.github.io/likelihood.model/articles/getting-started.html)
  goes deeper on the likelihood calculus generics.
- [`vignette("algebraic-mle-integration", package = "likelihood.model")`](https://queelius.github.io/likelihood.model/articles/algebraic-mle-integration.html)
  covers the `fisher_mle` class and how it inherits from `mle`.
- [`vignette("strategy-design", package = "compositional.mle")`](https://queelius.github.io/compositional.mle/articles/strategy-design.html)
  details the composition operators and why SICP-style closure matters
  for solvers.
