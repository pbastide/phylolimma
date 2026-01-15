# Transform data matrix

Multiply by inverse cholesky to whiten the data

## Usage

``` r
transform_data_tree(C_tree, y_data)
```

## Arguments

- C_tree:

  Cholesky of the tree structure obtained through
  [`get_chol_tree`](https://pbastide.github.io/phyloDE/reference/get_chol_tree.md)

- y_data:

  data matrix

## Value

transformed data matrix
