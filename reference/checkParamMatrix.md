# Check Matrix Parameter

Check that the parameters are compatible with the tree. Throws an error
if not.

## Usage

``` r
checkParamMatrix(x, name, tree, transpose = FALSE)
```

## Arguments

- x:

  matrix of parameters being tested.

- name:

  name of the parameter.

- tree:

  A phylogenetic tree with n tips.

- transpose:

  Should the transpose of the matrix be taken ? Default to FALSE.
