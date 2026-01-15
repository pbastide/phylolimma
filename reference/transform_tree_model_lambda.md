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
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html). It must be
  either a tree with tips having the same names as the columns of
  `object` (including replicates), or a tree such that tip labels match
  with species names in \`col_species\`.

- measurement_error:

  a logical value indicating whether there is measurement error. Default
  to `TRUE`. See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

## Value

The transformed tree.
