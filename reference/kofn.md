# Create a k-out-of-n system estimation model

Constructs a likelihood model for component lifetime estimation from
k-out-of-n system data. The system fails when `k` components have
failed: `k=1` is a series system (one failure kills it), `k=m` is a
parallel system (all must fail).

## Usage

``` r
kofn(
  k = 1L,
  m = 2L,
  component = dfr_exponential(),
  method = "mle",
  lifetime = "t",
  omega = "omega",
  lifetime_upper = "t_upper"
)
```

## Arguments

- k:

  System parameter: system fails when k components have failed. `k=1` is
  series, `k=m` is parallel.

- m:

  Number of components.

- component:

  A `dfr_dist` prototype from flexhaz specifying the shared component
  distribution family. Currently supported:
  [`dfr_exponential`](https://queelius.github.io/flexhaz/reference/dfr_exponential.html)
  and
  [`dfr_weibull`](https://queelius.github.io/flexhaz/reference/dfr_weibull.html).
  Parameter values on the prototype are ignored; they will be estimated
  from data.

- method:

  Estimation method: `"mle"` (direct MLE) or `"em"` (EM algorithm,
  Weibull only).

- lifetime:

  Column name for system lifetime (default `"t"`).

- omega:

  Column name for observation type (default `"omega"`).

- lifetime_upper:

  Column name for interval upper bound (default `"t_upper"`).

## Value

An S3 object of class
`c("exp_kofn"/"wei_kofn", "kofn", "likelihood_model")`.

## Details

This model satisfies the `likelihood_model` concept from the
`likelihood.model` package by providing methods for
[`loglik`](https://queelius.github.io/likelihood.model/reference/loglik.html),
[`score`](https://queelius.github.io/likelihood.model/reference/score.html),
and
[`hess_loglik`](https://queelius.github.io/likelihood.model/reference/hess_loglik.html).

The class hierarchy is:

- `"exp_kofn"` or `"wei_kofn"` (distribution-specific dispatch)

- `"kofn"` (shared methods)

- `"likelihood_model"` (generic inference infrastructure)

The concrete subclass is determined by the type of `component`. The
current closed-form machinery is homogeneous: all m components share the
same distribution family, and `component` acts as a prototype for that
family.

For non-k-of-n topologies (bridges, arbitrary coherent systems), use
[`coherent_dist`](https://rdrr.io/pkg/dist.structure/man/coherent_dist.html)
or one of the topology shortcuts in `dist.structure` directly. kofn is
exclusively for the k-out-of-n family.

## Examples

``` r
# Parallel system with 3 exponential components (k = m)
model <- kofn(k = 3, m = 3, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
print(model)
#> Error: object 'model' not found

# Series system (k = 1)
model_series <- kofn(k = 1, m = 4, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"

# Weibull parallel system with EM estimation
model_wei <- kofn(k = 2, m = 2, component = dfr_weibull(), method = "em")
#> Error in dfr_weibull(): could not find function "dfr_weibull"
```
