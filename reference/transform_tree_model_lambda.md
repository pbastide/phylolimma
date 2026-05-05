# Get lambda transformed tree

Compute the transformed tree using
[`transf.branch.lengths`](https://rdrr.io/pkg/phylolm/man/transf.branch.lengths.html).

## Usage

``` r
transform_tree_model_lambda(phy, phyfit, measurement_error)
```

## Arguments

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

- measurement_error:

  a logical value indicating whether there is measurement error, or
  individual independent (non phylogenetic) variation among samples.
  Default to `TRUE`. Setting this to `FALSE` can give unexpected
  results, except for the "lambda" model. See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

## Value

The transformed tree.
