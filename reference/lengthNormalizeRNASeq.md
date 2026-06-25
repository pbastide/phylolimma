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

  the RNASeq count matrix. Rows and columns should be named.

- lengthMatrix:

  the associated length matrix. Should have the same dimensions as
  `countMatrix`, with the same names.

- normalisationFactor:

  normalization factors to scale the raw library sizes, as computed e.g.
  by
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

- lengthNormalization:

  one of "none", "TPM" (default) or "RPKM". See details.

- dataTransformation:

  one of "log2" (default), "asin(sqrt)" or "sqrt." See details.

## Value

A matrix of normalized and transformed counts, with the same dimensions
as `countMatrix`.

## Details

The `lengthMatrix` is used to normalize the counts, using one of the
following formulas:

- `lengthNormalization="none"` : \\CPM\_{gi} = \frac{N\_{gi} + 0.5}{NF_i
  \times \sum\_{g} N\_{gi} + 1} \times 10^6\\

- `lengthNormalization="TPM"` : \\TPM\_{gi} = \frac{(N\_{gi} + 0.5) /
  L\_{gi}}{NF_i \times \sum\_{g} N\_{gi}/L\_{gi} + 1} \times 10^6\\

- `lengthNormalization="RPKM"` : \\RPKM\_{gi} = \frac{(N\_{gi} + 0.5) /
  L\_{gi}}{NF_i \times \sum\_{g} N\_{gi} + 1} \times 10^9\\

where \\N\_{gi}\\ is the count for gene g and sample i, \\L\_{gi}\\ is
the length of gene g in sample i, and \\NF_i\\ is the normalization for
sample i stored in `normalisationFactor`.

The function specified by the `dataTransformation` is then applied to
the normalized count matrix.

The "\\+0.5\\" is taken from Law et al 2014, and dropped from the
normalization when the transformation is something else than `log2`.

The "\\\times 10^6\\" and "\\\times 10^9\\" factors are omitted when the
`asin(sqrt)` transformation is taken, as \\asin\\ can only be applied to
real numbers smaller than 1.

## References

Law, C. W., Chen, Y., Shi, W. and Smyth, G. K. (2014), 'voom: precision
weights unlock linear model analysis tools for RNA-seq read counts',
Genome Biology 15(2), R29.

Bastide, P., Soneson, C., Stern, D. B., Lespinet, O. and Gallopin, M.
(2023), 'A Phylogenetic Framework to Simulate Synthetic Interspecies
RNA-Seq Data', Molecular Biology and Evolution 40(1), msac269.
