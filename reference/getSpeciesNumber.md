# Get the number of species

Compute the number of different species on a tree that possibly has
replicates coded as tips with zero length branches.

## Usage

``` r
getSpeciesNumber(phy, tol = .Machine$double.eps^(1/2))
```

## Arguments

- phy:

  a phylogentic tree, with possible replicates coded as tips with zero
  length branches.

- tol:

  a numeric value giving the tolerance to consider a branch length
  significantly greater than zero.

## Value

the number of different species in the tree
