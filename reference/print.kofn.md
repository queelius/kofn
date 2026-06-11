# Print method for kofn models

Displays a human-readable summary of the k-out-of-n system model,
including system type, component distribution, estimation method, and
column name conventions.

## Usage

``` r
# S3 method for class 'kofn'
print(x, ...)
```

## Arguments

- x:

  A `kofn` model object.

- ...:

  Additional arguments (ignored).

## Value

The model object, invisibly.

## Examples

``` r
print(kofn(k = 3, m = 3, component = dfr_exponential()))
#> Error in dfr_exponential(): could not find function "dfr_exponential"
```
