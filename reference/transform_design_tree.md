# Transform design matrix

Multiply by inverse cholesky to whiten the data

## Usage

``` r
transform_design_tree(C_tree, design)
```

## Arguments

- C_tree:

  Cholesky of the tree structure obtained through
  [`get_chol_tree`](https://pbastide.github.io/phyloDE/reference/get_chol_tree.md)

- design:

  design matrix

## Value

transformed design matrix
