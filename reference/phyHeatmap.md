# Heatmap with Phylogeny Structured Columns

This function uses the [`heatmap`](https://rdrr.io/r/stats/heatmap.html)
function to plot the (normalized) data, using the phylogenetic tree as
the column dendogram.

## Usage

``` r
phyHeatmap(
  object,
  design = NULL,
  coef = NULL,
  phy,
  scale = "none",
  ColSideColors = NULL,
  add.expr = NULL,
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

- coef:

  column number or column name specifying which coefficient or contrast
  of the linear model in the `design` matrix is of interest. If left
  `NULL`, defaults to the last column of `design`.

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

- scale:

  character indicating if the values on the heatmap should be centered
  and scaled in either the row direction or the column direction, or
  none. The default is "none", as the data is assumed to be already
  normalized. This scaling only affects the colour scale; it does not
  scale the original data. See documentation of
  [`heatmap`](https://rdrr.io/r/stats/heatmap.html).

- ColSideColors:

  (optional) character vector of length `ncol(object)` containing the
  color names for a horizontal side bar that may be used to annotate the
  columns of the heatmap. If left `NULL` and `design` is specified, the
  column of the `design` matrix corresponding to `coef` will be used.

- add.expr:

  expression that will be evaluated after the call to image. Can be used
  to add components to the plot. See documentation of
  [`heatmap`](https://rdrr.io/r/stats/heatmap.html). If left `NULL` and
  `design` is specified, the column of the `design` matrix corresponding
  to `coef` will be used to draw vertical lines.

- ...:

  further arguments to be passed to
  [`heatmap`](https://rdrr.io/r/stats/heatmap.html). Argument `Colv`
  defaults to the phylogeny, and cannot be overwritten.

## Value

Invisibly, a list, see
[`heatmap`](https://rdrr.io/r/stats/heatmap.html).
