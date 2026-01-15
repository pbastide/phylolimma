# Normalize RNASeq count data using gene lengths

This function normalizes a count matrix, using the matrix of length, and
an appropriate transformation.

## Usage

``` r
lengthNormalizeRNASeq(
  countMatrix,
  lengthMatrix = NULL,
  normalisationFactor = NULL,
  lengthNormalization = c("TPM", "RPKM", "none"),
  dataTransformation = c("log2", "sqrt", "asin(sqrt)")
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

- lengthNormalization:

  one of "none" (no correction), "TPM" (default) or "RPKM". See details.

- dataTransformation:

  one of "log2", "asin(sqrt)" or "sqrt." See details.

## Value

A matrix of normalized and transformed counts, with the same dimensions
as `countMatrix`.

## Details

The normalization procedures are:

- `none`::

  No length normalization.

- `TPM`::

  TODO

- `RPKM`::

  TODO
