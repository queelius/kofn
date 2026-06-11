# Assumptions for exponential k-out-of-n model

Returns a character vector listing the assumptions made by the
exponential k-out-of-n likelihood model.

## Usage

``` r
# S3 method for class 'exp_kofn'
assumptions(model, ...)
```

## Arguments

- model:

  An `exp_kofn` object created by
  [`kofn()`](https://queelius.github.io/kofn/reference/kofn.md).

- ...:

  Additional arguments (ignored).

## Value

Character vector of assumptions.

## Examples

``` r
model <- kofn(k = 3, m = 3, component = dfr_exponential())
#> Error in dfr_exponential(): could not find function "dfr_exponential"
assumptions(model)
#> Error: object 'model' not found
```
