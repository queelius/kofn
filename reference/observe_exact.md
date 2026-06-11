# Exact observation scheme (no censoring)

Creates an observation functor that records the exact failure time with
no censoring. This is the trivial case — included for completeness and
for explicit use in Fisher information comparisons.

## Usage

``` r
observe_exact()
```

## Value

Observation functor: `function(t_true)` returning a list with `t` (exact
time), `omega` (`"exact"`), and `t_upper` (`NA`).

## Examples

``` r
obs <- observe_exact()
obs(42.5)  # list(t = 42.5, omega = "exact", t_upper = NA)
#> $t
#> [1] 42.5
#> 
#> $omega
#> [1] "exact"
#> 
#> $t_upper
#> [1] NA
#> 
```
