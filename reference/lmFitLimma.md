# Fit using limma

Fit using limma

## Usage

``` r
lmFitLimma(y_trans, design_trans, ...)
```

## Arguments

- y_trans:

  A matrix data object containing normalized and phylogeny transformed
  expression values, with rows corresponding to genes and columns to
  samples (species).

- design_trans:

  the phylogeny transformed design matrix of the experiment, with rows
  corresponding to samples and columns to coefficients to be estimated.
  Defaults to the unit vector (intercept).

- ...:

  further parameters to be passed to
  [`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html).

## Value

A list with all the results.
