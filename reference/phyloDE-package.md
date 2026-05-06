# The `phyloDE` Package for Differential Expression Analysis

The `phyloDE` package fits linear model on inter-species gene expression
data, combining Phylogenetic Comparative Methods implemented in
[`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) with moderated
statistics tailored for gene expression implemented in
[`limma`](https://rdrr.io/pkg/limma/man/01Introduction.html).

With a design matrix that expresses a grouping conditions at the tip of
a phylogeny, `phyloDE` can perform Differential Expression analysis.

The main function of the package is
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md),
that inherits from the interfaces of both
[`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) and
[`limma`](https://rdrr.io/pkg/limma/man/01Introduction.html).

## See also

Useful links:

- <https://pbastide.github.io/phyloDE/>

## Author

Paul Bastide, Mélina Gallopin, Arnaud Liehrmann
