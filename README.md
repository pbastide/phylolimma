phyloDE
===============

<!-- badges: start -->
[![R-CMD-check](https://github.com/pbastide/phyloDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pbastide/phyloDE/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/pbastide/phyloDE/graph/badge.svg)](https://app.codecov.io/gh/pbastide/phyloDE)
<!-- badges: end -->

Functions to test for differential expression between conditions at the tip of a phylogeny, 
combining Phylogenetic Comparative Methods implemented in 
[`phylolm`](https://CRAN.R-project.org/package=phylolm)
with moderated statistics tailored for gene expression implemented in
[`limma`](https://bioconductor.org/packages/limma/).

## Installation

The development version can be installed from GitHub using the `remotes` package:
```R
install.packages("remotes")
remotes::install_github(repo = "pbastide/phyloDE")
```

## Documentation

See package documentation (references and vignettes) here: https://pbastide.github.io/phyloDE/.

(Built with [`pkgdown`](https://github.com/r-lib/pkgdown)).
