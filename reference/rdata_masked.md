# Generate masked k-out-of-n data

Returns a closure that generates system lifetime data with candidate
failed sets. The true failed set (the k components that failed by T_sys)
is always included in the candidate set. Additional non-failed
components are included independently with probability `p_mask`.

## Usage

``` r
rdata_masked(model, candset = "c", ...)
```

## Arguments

- model:

  A `kofn` model object.

- candset:

  Prefix for candidate set columns (default `"c"`).

- ...:

  Additional arguments (currently unused).

## Value

A function `function(theta, n, p_mask = 0, observe = NULL)` returning a
data frame with columns for system lifetime, observation type, and
Boolean candidate set indicators.

## Examples

``` r
model <- kofn(k = 2, m = 4, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
gen <- rdata_masked(model)
#> Error: object 'model' not found
set.seed(42)
df <- gen(theta = c(1, 0.8, 0.6, 0.4), n = 10, p_mask = 0.3)
#> Error in gen(theta = c(1, 0.8, 0.6, 0.4), n = 10, p_mask = 0.3): could not find function "gen"
head(df)
#>                                               
#> 1 function (x, df1, df2, ncp, log = FALSE)    
#> 2 {                                           
#> 3     if (missing(ncp))                       
#> 4         .Call(C_df, x, df1, df2, log)       
#> 5     else .Call(C_dnf, x, df1, df2, ncp, log)
#> 6 }                                           
```
