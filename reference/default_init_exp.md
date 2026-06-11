# Compute default initial values for exponential rates

Uses the method-of-moments estimator based on the harmonic number
relationship for i.i.d. parallel exponentials:
E[T](https://rdrr.io/r/base/logical.html) = H_m / lambda. Returns a
spread of initial values to break permutation symmetry.

## Usage

``` r
default_init_exp(df, m, lifetime, omega, lifetime_upper)
```

## Arguments

- df:

  Data frame with lifetime observations.

- m:

  Number of components.

- lifetime:

  Column name for lifetime.

- omega:

  Column name for observation type.

- lifetime_upper:

  Column name for interval upper bound.

## Value

Numeric vector of length `m` with initial rate estimates.
