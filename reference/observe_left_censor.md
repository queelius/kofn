# Left-censoring observation scheme

Creates an observation functor that applies left-censoring at time
`tau`. Systems that fail after `tau` are observed exactly; systems that
fail before `tau` are left-censored (we only know T \< tau).

## Usage

``` r
observe_left_censor(tau)
```

## Arguments

- tau:

  Censoring time (positive numeric).

## Value

Observation functor: `function(t_true)` returning a list with `t`,
`omega` (`"left"` or `"exact"`), and `t_upper` (`NA`).

## Examples

``` r
obs <- observe_left_censor(tau = 10)
obs(5)   # left:  list(t = 10, omega = "left", t_upper = NA)
#> $t
#> [1] 10
#> 
#> $omega
#> [1] "left"
#> 
#> $t_upper
#> [1] NA
#> 
obs(15)  # exact: list(t = 15, omega = "exact", t_upper = NA)
#> $t
#> [1] 15
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
#> 
```
