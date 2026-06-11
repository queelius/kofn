# Periodic inspection observation scheme

Creates an observation functor for periodic inspections at intervals of
width `delta`. The failure time is known only to lie within the
inspection interval containing it (interval-censored). Systems surviving
past `tau` are right-censored.

## Usage

``` r
observe_periodic(delta, tau = Inf)
```

## Arguments

- delta:

  Inspection interval width (positive numeric).

- tau:

  Study end time / right-censoring time (positive numeric or `Inf` for
  no right-censoring).

## Value

Observation functor: `function(t_true)` returning a list with `t`
(interval lower bound or `tau`), `omega` (`"interval"` or `"right"`),
and `t_upper` (interval upper bound or `NA`).

## Examples

``` r
obs <- observe_periodic(delta = 10, tau = 100)
obs(25)   # interval: list(t = 20, omega = "interval", t_upper = 30)
#> $t
#> [1] 20
#> 
#> $omega
#> [1] "interval"
#> 
#> $t_upper
#> [1] 30
#> 
obs(150)  # right:    list(t = 100, omega = "right", t_upper = NA)
#> $t
#> [1] 100
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA
#> 
```
