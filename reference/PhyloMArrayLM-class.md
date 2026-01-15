# Class PhyloMArrayLM

Extension of class
[`MArrayLM-class`](https://rdrr.io/pkg/limma/man/marraylm.html) from
package `limma` for phylogenetically correlated regressions of gene-wise
linear models.

It is a simple list-based S4 class. Objects are normally created by
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).
Additional components are added by
[`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html).

## Components

`PhyloMArrayLM` objects do not contain any slots (apart from .Data) but
they should contain the same list components than
[`MArrayLM-class`](https://rdrr.io/pkg/limma/man/marraylm.html).

In addition, they contain the following tree-specific components:

- `phy`::

  The phylogenetic tree used for the regression, of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- `modelphy`::

  The phylogenetic model of trait evolution, argument call in
  [`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).

- `measurement_error`::

  Boolean, TRUE if there is additional measurement error, argument call
  in
  [`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).

- `phy_trans`::

  List of the transformed phylogenetic trees obtained from a fit using
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) on each
  gene.

- `optpar`::

  Vector of the `optpar` obtained from a fit using
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) on each
  gene.

- `lambda_error`::

  Vector of the `lambda_error` parameters obtained from a fit using
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) on each
  gene.

- `sigma2_phy`::

  Vector of the phylogenetic `sigma2` parameters obtained from a fit
  using [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) on
  each gene.

- `sigma2_error`::

  Vector of the `sigma2_error` parameters obtained from a fit using
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) on each
  gene.

## See also

[`MArrayLM-class`](https://rdrr.io/pkg/limma/man/marraylm.html),
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).
