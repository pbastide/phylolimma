# Class ConsensusTreeModel

A simple list-based S4 class, that contains a tree and associated
parameters, obtained through function
[`phylogeneticCorrelations`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md).

## Components

`ConsensusTreeModel` objects do not contain any slots (apart from .Data)
but they should contain the following list components:

\#'

- `tree` the transformed consensus tree

- `params` the associated consensus parameters.

&nbsp;

- `tree`::

  the transformed consensus tree used to define the correlation
  structure.

- `params`::

  the parameters associated with the tree, including regularized values
  of the parameters.

## See also

[`phylogeneticCorrelations`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md),
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md),
[`PhyloMArrayLM-class`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md).
