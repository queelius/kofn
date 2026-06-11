# Right-censoring observation scheme

Creates an observation functor that applies right-censoring at time
`tau`. Systems that fail before `tau` are observed exactly; systems
surviving past `tau` are right-censored.

## Usage

``` r
observe_right_censor(tau = Inf)
```

## Arguments

- tau:

  Censoring time (positive numeric, or `Inf` for no censoring).

## Value

Observation functor: `function(t_true)` returning a list with `t`,
`omega` (`"exact"` or `"right"`), and `t_upper` (`NA`).

## Examples

``` r
obs <- observe_right_censor(tau = 100)
obs(50)   # exact:  list(t = 50, omega = "exact", t_upper = NA)
#> $t
#> [1] 50
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
#> 
obs(150)  # right:  list(t = 100, omega = "right", t_upper = NA)
#> $t
#> [1] 100
#> 
#> $omega
#> [1] "right"
#> 
#> $t_upper
#> [1] NA
#> 

# No censoring
obs_all <- observe_right_censor(tau = Inf)
obs_all(999)  # exact: list(t = 999, omega = "exact", t_upper = NA)
#> $t
#> [1] 999
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
#> 
```
