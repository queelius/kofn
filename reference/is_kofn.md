# Test if an object is a kofn model

Test if an object is a kofn model

## Usage

``` r
is_kofn(x)
```

## Arguments

- x:

  An object to test.

## Value

Logical indicating whether `x` inherits from `"kofn"`.

## Examples

``` r
is_kofn(kofn(k = 3, m = 3, component = dfr_exponential()))
#> Error in dfr_exponential(): could not find function "dfr_exponential"
is_kofn(42)
#> [1] FALSE
```
