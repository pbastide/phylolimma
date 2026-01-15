# Trimmed mean

Trimmed mean, with asymetric trim up and down. Code adapted from
\`mean.default\`, with a different trim value for the lower and upper
ends.

## Usage

``` r
mean_trim(x, trim = 0.15, na.rm = FALSE, ...)
```

## Arguments

- x:

  a numeric vector.

- trim:

  a vector of size two, with the proportion of lower and upper values to
  trim. If a single value is provided, it will be used for lower and
  upper trim.

- na.rm:

  whether to remove \`NA\`s before the computation.

- ...:

  further arguments passed to or from other methods.

## Value

The trimmed arithmetic mean.
