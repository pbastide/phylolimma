# Package index

## Phylogenetic Differential Expression

Functions to test for differential expression using moderated
phylogenetic linear regression.

- [`phyloDE-package`](https://pbastide.github.io/phyloDE/reference/phyloDE-package.md)
  [`phyloDE`](https://pbastide.github.io/phyloDE/reference/phyloDE-package.md)
  :

  The `phyloDE` Package for Differential Expression Analysis

- [`phylolmFit()`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md)
  : Phylogenetic Linear Model for Gene Expression Analysis

- [`phylogeneticCorrelations()`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md)
  : Phylogenetic Correlation using a Consensus Tree

## S4 Classes

Classes containing the result of a fit and associated helper functions.

- [`PhyloMArrayLM-class`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md)
  : Class PhyloMArrayLM
- [`ConsensusTreeModel-class`](https://pbastide.github.io/phyloDE/reference/ConsensusTreeModel-class.md)
  : Class ConsensusTreeModel
- [`getParameters()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`plotParameters()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`show(`*`<PhyloMArrayLM>`*`)`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`show(`*`<ConsensusTreeModel>`*`)`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`logLikelihood()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`consensusTree()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  : Methods for classes PhyloMArrayLM and ConsensusTreeModel

## Helper Functions

Function to pre-process the tree or the data.

- [`phyHeatmap()`](https://pbastide.github.io/phyloDE/reference/phyHeatmap.md)
  : Heatmap with Phylogeny Structured Columns
- [`addReplicatesOnTree()`](https://pbastide.github.io/phyloDE/reference/addReplicatesOnTree.md)
  : Add replicates to a tree
- [`lengthNormalizeRNASeq()`](https://pbastide.github.io/phyloDE/reference/lengthNormalizeRNASeq.md)
  : Normalize RNASeq count data using gene lengths
- [`rhoFromAlpha()`](https://pbastide.github.io/phyloDE/reference/rhoFromAlpha.md)
  [`alphaFromRho()`](https://pbastide.github.io/phyloDE/reference/rhoFromAlpha.md)
  : Compute rho parameter
- [`getBoundsSelectionStrength()`](https://pbastide.github.io/phyloDE/reference/getBoundsSelectionStrength.md)
  : Get bounds on alpha for an OU
- [`getMinError()`](https://pbastide.github.io/phyloDE/reference/getMinError.md)
  : Get Lower Bound on sigma2_error

## Dataset

Crayfish dataset

- [`crayfish`](https://pbastide.github.io/phyloDE/reference/crayfish.md)
  : Crayfish RNA-Seq dataset
