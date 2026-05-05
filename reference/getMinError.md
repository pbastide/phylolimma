# Get Lower Bound on sigma2_error

Find reasonable lower bound on the `sigma2_error` parameter when fitted
with [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html).

## Usage

``` r
getMinError(phy, tol = (.Machine$double.eps)^0.5)
```

## Arguments

- phy:

  a phylogenetic tree.

- tol:

  the numerical tolerance.

## Value

The minimum value for `sigma2_error`.

## Details

The minimum value must be high enough so that it can be numerically
distinguished from zero. This is important for trees with several
samples per species, as the limit case `sigma2_error = 0` is degenerate,
as it means that all the samples in one species have the exact same
value.

Default to \\tol \* h\\, where \\h\\ is the total height of the tree.
Note that if an OU is used in the fit, then this value is updated to
match the transformed tree height.
