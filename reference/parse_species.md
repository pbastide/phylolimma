# Add replicates to a tree

Utility function to add replicates to a tree, as tips with zero length
branches.

## Usage

``` r
parse_species(tree, ids, pattern = "(_|\\.).*$")
```

## Arguments

- tree:

  A phylogenetic tree with n tips.

- ids:

  a vector of sample ids.

- pattern:

  a regular expression to find species from sample names. Default to
  removing everything after a dot or underscore.

## Value

A vector of the same length as \`ids\`, with the species of each sample.
