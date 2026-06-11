# Interval-censoring observation scheme

Creates an observation functor that places all failure times into a
fixed window `[a, b)`. Systems failing within the window are
interval-censored; systems failing outside are observed exactly.

## Usage

``` r
observe_interval_censor(a, b)
```

## Arguments

- a:

  Lower bound of censoring window (non-negative numeric).

- b:

  Upper bound of censoring window (positive numeric, `b > a`).

## Value

Observation functor: `function(t_true)` returning a list with `t`,
`omega` (`"interval"` or `"exact"`), and `t_upper`.

## Details

For regular grids, use
[`observe_periodic`](https://queelius.github.io/kofn/reference/observe_periodic.md)
instead.

## Examples

``` r
obs <- observe_interval_censor(a = 5, b = 10)
obs(7)   # interval: list(t = 5, omega = "interval", t_upper = 10)
#> $t
#> [1] 5
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 10
#> 
obs(3)   # exact:    list(t = 3, omega = "exact", t_upper = NA)
#> $t
#> [1] 3
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
#> 
obs(12)  # exact:    list(t = 12, omega = "exact", t_upper = NA)
#> $t
#> [1] 12
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
#> 
```
