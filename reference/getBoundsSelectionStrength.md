# Get Bounds on alpha

Find reasonable bounds on the `alpha` parameter.

## Usage

``` r
getBoundsSelectionStrength(
  phy,
  relative_half_life_min = 1e-04,
  relative_half_life_max = 10000
)
```

## Arguments

- phy:

  a phylogenetic tree.

- relative_half_life_min:

  optional minimal half life relative to tree height

- relative_half_life_max:

  optional maximal half life relative to tree height

## Value

A vector with lower and upper bounds for alpha

## Details

This functions tries to find reasonable bounds on the \\\alpha\\
parameter of an OU process by using the scaled phylogenetic half-life
\\t\_{1/2} = \log(2) / \alpha / h\\, where \\h\\ is the total height of
the tree. If \\t\_{1/2} = D)\\, it means that the trait will need a time
\\D \times h\\ to cover half the distance to the optimal value (see
Hansen, 1997). Small values of \\D\\ means high selection pressure
(large \\\alpha\\), while large values of \\D\\ means low selection
pressure (small \\\alpha\\).

The default maximum value for \\D\\ is `relative_half_life_max = 10000`
(selection is week and the process looks like a BM).

The default minimum value for \\D\\ is `relative_half_life_min = 0.0001`
(selection is strong and tips are only weakly correlated).

The function makes sure that the maximum \\\alpha\\ value associated
with `relative_half_life_min` does not lead to underflow errors when
rescaling the tree.

## References

Hansen, T. F. (1997). Stabilizing Selection and the Comparative Analysis
of Adaptation. Evolution, 51(5) :1341.
