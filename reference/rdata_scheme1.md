# Generate Scheme 1 (periodic inspection) data

Returns a closure that generates k-out-of-n system data under periodic
inspection. Each observation yields the exact system failure time and
interval-censored component failure times.

## Usage

``` r
rdata_scheme1(model, ...)
```

## Arguments

- model:

  A `kofn` model object (exponential or Weibull).

- ...:

  Additional arguments (currently unused).

## Value

A function `function(theta, n, delta)` returning a data frame with
columns:

- `t`:

  System failure time (exact).

- `comp_lower_j`:

  Lower bound of inspection interval for component j.

- `comp_upper_j`:

  Upper bound of inspection interval for component j.

The data frame has attributes `comp_times` (true component times),
`delta` (inspection interval), and `par` (true parameters).

## Details

Note: The Scheme 1 likelihood uses a composite approximation that treats
the system density and component interval contributions as independent.
This works well for k \>= 2 but is unreliable for series systems (k =
1), where surviving components' intervals should be conditioned on
survival past T_sys. For series systems, use the `maskedcauses` package
instead.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
gen <- rdata_scheme1(model)
#> Error: object 'model' not found
set.seed(1)
df <- gen(theta = c(1, 2), n = 20, delta = 1.0)
#> Error in gen(theta = c(1, 2), n = 20, delta = 1): could not find function "gen"
head(df)
#>                                               
#> 1 function (x, df1, df2, ncp, log = FALSE)    
#> 2 {                                           
#> 3     if (missing(ncp))                       
#> 4         .Call(C_df, x, df1, df2, log)       
#> 5     else .Call(C_dnf, x, df1, df2, ncp, log)
#> 6 }                                           
```
