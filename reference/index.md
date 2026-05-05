# Package index

## Phylogenetic Differential Expression

Functions to test for differential expression using moderated
phylogenetic linear regression.

- [`phyloDE-package`](https://pbastide.github.io/phyloDE/reference/phyloDE-package.md)
  [`phyloDE`](https://pbastide.github.io/phyloDE/reference/phyloDE-package.md)
  : TODO
- [`phylolmFit()`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md)
  : Phylogenetic Linear Model for Gene Expression Analysis
- [`phylogeneticCorrelations()`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md)
  : Phylogenetic correlation using a consensus tree

## PhyloMArrayLM Class

Class containing the result of a fit and associated helper functions.

- [`PhyloMArrayLM-class`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md)
  : Class PhyloMArrayLM
- [`getParameters()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`plotParameters()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`show(`*`<PhyloMArrayLM>`*`)`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`logLikelihood()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  [`consensusTree()`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLMMethods.md)
  : Methods for class PhyloMArrayLM

## Helper Functions

Function to pre-process the tree or the data.

- [`addReplicatesOnTree()`](https://pbastide.github.io/phyloDE/reference/addReplicatesOnTree.md)
  : Add replicates to a tree
- [`lengthNormalizeRNASeq()`](https://pbastide.github.io/phyloDE/reference/lengthNormalizeRNASeq.md)
  : Normalize RNASeq count data using gene lengths
- [`rhoFromAlpha()`](https://pbastide.github.io/phyloDE/reference/rhoFromAlpha.md)
  [`alphaFromRho()`](https://pbastide.github.io/phyloDE/reference/rhoFromAlpha.md)
  : Compute rho parameter

## Dataset

Crayfish dataset

- [`crayfish`](https://pbastide.github.io/phyloDE/reference/crayfish.md)
  : Crayfish RNA-Seq dataset

## Not Implemented

Functions that apply to MArrayLM-class but not implemented for
PhyloMArrayLM-class

- [`treat()`](https://pbastide.github.io/phyloDE/reference/treat.md)
  [`decideTests()`](https://pbastide.github.io/phyloDE/reference/treat.md)
  [`classifyTestsF()`](https://pbastide.github.io/phyloDE/reference/treat.md)
  : Not Implemented Functions
