# Determinant of a positive-definite Hessian, or NA

Computes the numerical Hessian, checks finite-ness and positive
definiteness, and returns `det(H)` or `NA`.

## Usage

``` r
safe_hessian_det(neg_ll, par)
```

## Arguments

- neg_ll:

  Negative log-likelihood function.

- par:

  Parameter vector.

## Value

Scalar determinant, or `NA_real_`.
