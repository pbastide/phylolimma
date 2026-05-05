# Phylogenetic correlation using a consensus tree

This function applies
[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) to the normalized
data, in order to take the phylogeny into account. TODO: explain more.

## Usage

``` r
phylogeneticCorrelations(
  object,
  design = NULL,
  phy,
  col_species = NULL,
  model = c("OUfixedRoot", "BM", "lambda"),
  measurement_error = TRUE,
  trim = c(0.25, 0.05),
  REML = TRUE,
  ncores = 1,
  ...
)
```

## Arguments

- object:

  a matrix data object containing normalized expression values, with
  rows corresponding to genes and columns to samples (species).

- design:

  the design matrix of the experiment, with rows corresponding to
  samples and columns to coefficients to be estimated. Defaults to the
  unit vector (intercept).

- phy:

  an object of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html), representing
  the phylogenetic relationships between the species. It must be dated
  and ultrametric. If the column names of `object` follow the pattern
  `SpeciesName_SampleId` or `SpeciesName.SampleId`, an automatic
  matching of the samples on the tip of the tree is performed.
  Otherwise, the tree tip labels must match with species names in
  `col_species` (see below). The tip labels of the tree can also match
  exactly the names as the columns of `object`, so that the tree
  directly includes all the replicates.

- col_species:

  a character vector with same length as there are columns in the
  expression matrix, specifying the species for the corresponding
  column. If left `NULL` (the default), an automatic parsing of species
  names with sample ids is attempted.

- model:

  the phylogenetic model used to correct for the phylogeny. Must be one
  of "OUfixedRoot" (the default), "BM", or "lambda". See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

- measurement_error:

  a logical value indicating whether there is measurement error, or
  individual independent (non phylogenetic) variation among samples.
  Default to `TRUE`. Setting this to `FALSE` can give unexpected
  results, except for the "lambda" model. See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

- trim:

  a vector of size two, with the fraction of observations to be trimmed
  from the lower and upper ends of `atanh(all.lambdas)` and `atanh(rho)`
  when computing the trimmed mean. If a single value is provided, it is
  recycled as a vector of size two. Default to `c(0.25, 0.05)`. See also
  the `trim` argument in
  [`duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html).

- REML:

  Use REML (default) or ML for estimating the parameters.

- ncores:

  number of cores to use for parallel computation. Default to 1 (no
  parallel computation).

- ...:

  further parameters to be passed to
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html).

## Value

An object of class `TransTree-class`, with list components:

- `tree` the transformed consensus tree

- `params` the associated consensus parameters.

This consensus tree defines a correlation structure, and can be passed
on to [`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html).

## See also

[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html),
[`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html),
[`duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html)

[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html),
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md),
[`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html),
[`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html)
