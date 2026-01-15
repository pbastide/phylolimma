# Check the tree

Check the tree

## Usage

``` r
check_tree(phy, y, col_species)
```

## Arguments

- phy:

  an object of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html). It must be
  either a tree with tips having the same names as the columns of
  `object` (including replicates), or a tree such that tip labels match
  with species names in \`col_species\`.

- col_species:

  a character vector with same length as columns in the expression
  matrix, specifying the species for the corresponding column. If left
  \`NULL\`, an automatic parsing of species names with sample ids is
  attempted.

## Value

the correctly formatted tree
