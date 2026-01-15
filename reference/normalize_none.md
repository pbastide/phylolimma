# Normalize RNASeq count data

Apply standard CPM, with no length normalization.

## Usage

``` r
normalize_none(countMatrix, normalisationFactor, dataTransformation)
```

## Arguments

- countMatrix:

  The RNASeq count matrix. Rows and columns should be named.

- normalisationFactor:

  Normalization factors to scale the raw library sizes, as computed e.g.
  by
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

- dataTransformation:

  one of "log2", "asin(sqrt)" or "sqrt." See details.

## Value

A matrix of normalized count, with the same dimensions as `countMatrix`.
