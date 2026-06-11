# Assumptions for Weibull k-out-of-n system model

Returns a list of assumptions made by the model.

## Usage

``` r
# S3 method for class 'wei_kofn'
assumptions(model, ...)
```

## Arguments

- model:

  A `wei_kofn` object.

- ...:

  Additional arguments (currently unused).

## Value

A character vector of assumptions.

## Examples

``` r
model <- kofn(k = 2, m = 2, component = dfr_weibull())
#> Error in dfr_weibull(): could not find function "dfr_weibull"
assumptions(model)
#> Error: object 'model' not found
```
