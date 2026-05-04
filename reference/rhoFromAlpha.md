# Compute rho parameter

Function `rhoFromAlpha` computes the `rho` parameter, as defined in
Cornuault, 2023, Syst. Biol:
`rho = 1 - (1 - exp(-2 * alpha * t_tree)) / (2 * alpha * t_tree)` `rho`
can be interpreted as the percent decrease in trait variance caused by
the OU as compared to the variance expected under under BM. When alpha
goes to 0, `rho` goes to 0: the trait can be explained by a "neutral" BM
process. When alpha goes to Inf, `rho` goes to 1: the trait exhibits no
tree structure.

Function `alphaFromRho` computes the value of `alpha` from the value of
`rho`. This function can be expressed using the principal Lambert W
function. Here, it is computed by inverting function `rhoFromAlpha`
directly, using function
[`uniroot`](https://rdrr.io/r/stats/uniroot.html) within the
`alpha_bounds` bounds.

## Usage

``` r
rhoFromAlpha(alpha, t_tree, tol = .Machine$double.eps)

alphaFromRho(rho, t_tree, alpha_bounds = NULL)
```

## Arguments

- alpha:

  the selection strength of the process

- t_tree:

  the total height of the tree

- tol:

  tolerence value for calling `alpha = 0`.

- rho:

  the rho value

- alpha_bounds:

  lower and upper bounds on alpha values. If \`NULL\` (the default),
  default bounds are computed from the tree height.

## Value

Value of `rho` or `alpha`.
