# Generate random data from a Weibull k-out-of-n system

Returns a closure that generates system lifetime data with observation
scheme support. Component lifetimes are Weibull-distributed with
parameters specified by `theta` (interleaved shape/scale).

## Usage

``` r
# S3 method for class 'wei_kofn'
rdata(model, ...)
```

## Arguments

- model:

  A `wei_kofn` object.

- ...:

  Additional arguments (currently unused).

## Value

A function `function(theta, n, observe = NULL)` returning a data frame
with columns `t` (system lifetimes) and `omega` (observation type).
Attributes:

- `comp_times`:

  n x m matrix of component lifetimes.

- `par`:

  The parameter vector used for generation.

## Details

The interface matches
[`rdata.exp_kofn`](https://queelius.github.io/kofn/reference/rdata.exp_kofn.md):
the returned closure accepts an `observe` functor for observation
schemes (censoring, periodic inspection, etc.). The output includes an
`omega` column.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_weibull())
#> Error in dfr_weibull(): could not find function "dfr_weibull"
gen <- rdata(model)
#> Error: object 'model' not found
set.seed(42)
df <- gen(theta = c(1.5, 2.0, 2.0, 3.0), n = 50)
#> Error in gen(theta = c(1.5, 2, 2, 3), n = 50): could not find function "gen"
head(df)
#>                                               
#> 1 function (x, df1, df2, ncp, log = FALSE)    
#> 2 {                                           
#> 3     if (missing(ncp))                       
#> 4         .Call(C_df, x, df1, df2, log)       
#> 5     else .Call(C_dnf, x, df1, df2, ncp, log)
#> 6 }                                           

# With right-censoring
df_cens <- gen(c(1.5, 2.0, 2.0, 3.0), 50,
               observe = observe_right_censor(tau = 3))
#> Error in gen(c(1.5, 2, 2, 3), 50, observe = observe_right_censor(tau = 3)): could not find function "gen"
table(df_cens$omega)
#> Error: object 'df_cens' not found
```
