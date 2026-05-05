# Check the tree

Check the tree

## Usage

``` r
check_tree(phy, y, col_species)
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

- col_species:

  a character vector with same length as there are columns in the
  expression matrix, specifying the species for the corresponding
  column. If left `NULL` (the default), an automatic parsing of species
  names with sample ids is attempted.

## Value

the correctly formatted tree
