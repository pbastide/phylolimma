# Add replicates to a tree

Utility function to add replicates to a tree, as tips with zero length
branches.

## Usage

``` r
addReplicatesOnTree(
  tree,
  traits,
  species = "species",
  id = "id",
  eps = .Machine$double.eps^2
)
```

## Arguments

- tree:

  A phylogenetic tree with n tips.

- traits:

  A data frame containing at least two columns, one with sample ids, and
  one with species names for each samples.

- species:

  Name of the column containing species names. Default to "species".

- id:

  Name of the column containing samples ids. Default to "id".

- eps:

  A small number to add to terminal branch lengths to avoid true zeros.
  Default to `.Machine$double.eps^2`.

## Value

A phylogenetic tree with as many tips as the number of rows in `traits`,
and clusters of tips with zero branch lengths corresponding to
replicates.
