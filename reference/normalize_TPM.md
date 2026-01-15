# Normalize RNASeq count data

Apply TPM length normalization.

## Usage

``` r
normalize_TPM(
  countMatrix,
  lengthMatrix,
  normalisationFactor,
  dataTransformation
)
```

## Arguments

- countMatrix:

  The RNASeq count matrix. Rows and columns should be named.

- lengthMatrix:

  The associated length matrix. Should have the same dimensions as
  `countMatrix`, with the same names.

- normalisationFactor:

  Normalization factors to scale the raw library sizes, as computed e.g.
  by
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

- dataTransformation:

  one of "log2", "asin(sqrt)" or "sqrt." See details.

## Value

A matrix of normalized count, with the same dimensions as `countMatrix`.
