# Number of components in a kofn model

Returns the number of components `m` in the k-out-of-n system.

## Usage

``` r
# S3 method for class 'kofn'
ncomponents(x, ...)
```

## Arguments

- x:

  A `kofn` model object.

- ...:

  Additional arguments (ignored).

## Value

Integer number of components.

## Examples

``` r
ncomponents(kofn(k = 5, m = 5, component = dfr_exponential()))
#> Error in ncomponents(kofn(k = 5, m = 5, component = dfr_exponential())): could not find function "ncomponents"
```
