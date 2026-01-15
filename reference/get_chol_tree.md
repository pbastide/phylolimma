# Get Tree Normalizing Inverse Cholesky

Compute the whitening cholesky matrix.

## Usage

``` r
get_chol_tree(y_data, design, phy, model, measurement_error, REML, ncores, ...)
```

## Arguments

- y_data:

  A matrix data object containing normalized expression values, with
  rows corresponding to genes and columns to samples (species).

- design:

  the design matrix of the experiment, with rows corresponding to
  samples and columns to coefficients to be estimated. Defaults to the
  unit vector (intercept).

- phy:

  an object of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html). It must be
  either a tree with tips having the same names as the columns of
  `object` (including replicates), or a tree such that tip labels match
  with species names in \`col_species\`.

- model:

  the phylogenetic model used to correct for the phylogeny. Must be one
  of "BM", "lambda" or "OUfixedRoot". See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

- measurement_error:

  a logical value indicating whether there is measurement error. Default
  to `TRUE`. See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

- REML:

  Use REML (default) or ML for estimating the parameters.

- ncores:

  number of cores to use for parallel computation. Default to 1 (no
  parallel computation).

- ...:

  further parameters to be passed to
  [`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) or
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html).

## Value

The (list of) cholesky matrix of the tree structure.
