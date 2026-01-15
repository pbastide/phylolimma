# Get Min on sigma2_error

Find reasonable minimum on the `sigma2_error` parameter.

## Usage

``` r
getMinError(phy, tol = (.Machine$double.eps)^0.5)
```

## Arguments

- phy:

  a phylogenetic tree.

- tol:

  the numerical tolerance

## Value

The minimum value for sigma2_error

## Details

The minimum value must be high enough so that it can be numerically
distinguished from zero. Default to \\tol \* h\\, where \\h\\ is the
total height of the tree. If an OU is used, then this value is updated
to match the transformed tree height.
