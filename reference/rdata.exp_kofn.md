# Random data generation for exponential k-out-of-n model

Returns a closure that generates random system-level data from the
exponential k-out-of-n data-generating process.

## Usage

``` r
# S3 method for class 'exp_kofn'
rdata(model, ...)
```

## Arguments

- model:

  An `exp_kofn` object created by
  [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md).

- ...:

  Additional arguments (ignored).

## Value

A function `function(theta, n, observe = NULL)` returning a data frame
with columns `t` (system lifetime) and `omega` (observation type).
Latent component lifetimes, true critical component, and the true
parameters are stored as attributes.

## Details

Workflow:

1.  Generate i.i.d. exponential component lifetimes.

2.  Compute system lifetime as the (m - k + 1)-th order statistic.

3.  Apply the observation functor (exact observation by default).

## Examples

``` r
model <- kofn(k = 3, m = 3, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
gen <- rdata(model)
#> Error: object 'model' not found
set.seed(1)
df <- gen(theta = c(1, 2, 3), n = 20)
#> Error in gen(theta = c(1, 2, 3), n = 20): could not find function "gen"
head(df)
#>                                               
#> 1 function (x, df1, df2, ncp, log = FALSE)    
#> 2 {                                           
#> 3     if (missing(ncp))                       
#> 4         .Call(C_df, x, df1, df2, log)       
#> 5     else .Call(C_dnf, x, df1, df2, ncp, log)
#> 6 }                                           

# With right-censoring
df2 <- gen(c(1, 2, 3), n = 20, observe = observe_right_censor(tau = 2))
#> Error in gen(c(1, 2, 3), n = 20, observe = observe_right_censor(tau = 2)): could not find function "gen"
table(df2$omega)
#> Error: object 'df2' not found
```
