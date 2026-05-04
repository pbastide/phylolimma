# Phylogenetic Comparative Method using LIMMA

This function applies
[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) to the normalized
data, in order to take the phylogeny into account. TODO: explain more.

## Usage

``` r
phylolmFit(
  object,
  design = NULL,
  phy,
  col_species = NULL,
  model = c("OUfixedRoot", "BM", "lambda"),
  measurement_error = TRUE,
  use_consensus = TRUE,
  consensus_tree = NULL,
  REML = TRUE,
  ncores = 1,
  ...
)
```

## Arguments

- object:

  A matrix data object containing normalized expression values, with
  rows corresponding to genes and columns to samples (species).

- design:

  the design matrix of the experiment, with rows corresponding to
  samples and columns to coefficients to be estimated. Defaults to the
  unit vector (intercept).

- phy:

  an object of class
  [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html). It must be
  either a tree with tips having the same names as the columns of
  `object` (including replicates), or a tree such that tip labels match
  with species names in \`col_species\`.

- col_species:

  a character vector with same length as columns in the expression
  matrix, specifying the species for the corresponding column. If left
  \`NULL\`, an automatic parsing of species names with sample ids is
  attempted.

- model:

  the phylogenetic model used to correct for the phylogeny. Must be one
  of "BM", "lambda" or "OUfixedRoot". See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

- measurement_error:

  a logical value indicating whether there is measurement error. Default
  to `TRUE`. See
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html) for more
  details.

- use_consensus:

  If `TRUE`, one consensus tree is used to represent the correlation
  structure. see
  [`phylogeneticCorrelations`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md).
  If `FALSE`, each gene will use its own model parameters and will have
  its own correlation structure accordingly. Default to TRUE.

- consensus_tree:

  If not `NULL`, the consensus tree containing the correlation
  structure, result of function
  [`phylogeneticCorrelations`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md).
  If provided, arguments `phy`, `model` and `measurement_error` will be
  ignored.

- REML:

  Use REML (default) or ML for estimating the parameters.

- ncores:

  number of cores to use for parallel computation. Default to 1 (no
  parallel computation).

- ...:

  further parameters to be passed to
  [`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) or
  [`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html).

## Value

An object of class
[`PhyloMArrayLM-class`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md),
with list components `coefficients`, `stdev.unscaled`, `sigma` and
`df.residual`. These quantities take the phylogenetic model into
account. The object can be passed to
[`eBayes`](https://pbastide.github.io/phyloDE/reference/eBayes.md).

## Details

The default bounds on the phylogenetic parameters are the same as in
[`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html), except for
the `alpha` parameter of the OU.

## See also

[`lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html),
[`phylolm`](https://rdrr.io/pkg/phylolm/man/phylolm.html),
[`phylogeneticCorrelations`](https://pbastide.github.io/phyloDE/reference/phylogeneticCorrelations.md),
[`eBayes`](https://pbastide.github.io/phyloDE/reference/eBayes.md)

## Examples

``` r
## Simulate a tree
set.seed(1289)
ntips <- 20
tree <- ape::rphylo(ntips, 0.1, 0)
condition <- sample(c(0, 1), ntips, replace = TRUE)

## Simulate data with replicates
reps <- sample(1:5, ntips, replace = TRUE)
rep_ids <- make.unique(rep(tree$tip.label, times = reps), sep = "_")
ngenes <- 20
dat <- matrix(rnorm(sum(reps) * ngenes, 1, 0.5), nrow = ngenes)
rownames(dat) <- paste0("g", 1:ngenes)
colnames(dat) <- rep_ids

## Add differentially expressed genes
condition_reps <- rep(condition, times = reps)
ndiff <- 5
dat[1:ndiff, ] <- dat[1:ndiff, ] + rexp(ndiff, 1/2) %*% t(condition_reps)

## Design matrix
design <- cbind(rep(1, ncol(dat)), condition_reps)
colnames(design) <- c("(Intercept)", "condition")
rownames(design) <- rep_ids

## linear model fit
pfit <- phylolmFit(dat, design = design, phy = tree)
#> Loading required package: ape

## eBayes correction
pfit <- eBayes(pfit, trend = TRUE)
limma::topTable(pfit, coef = 2)
#>          logFC   AveExpr          t       P.Value     adj.P.Val           B
#> g2  25.2948752 7.6462148 194.650854  0.000000e+00  0.000000e+00 1190.049817
#> g5   9.0026076 3.1593224  72.724947 1.262101e-293 1.262101e-292  661.316179
#> g4   7.7178819 2.6802085  62.294317 4.560948e-259 3.040632e-258  582.396074
#> g1   0.7354200 0.7938382   5.525255  4.975259e-08  2.487629e-07    6.885335
#> g3   0.6213833 0.6894336   4.654513  4.028371e-06  1.611348e-05    2.641743
#> g15  0.1925638 0.6160945   1.446090  1.486921e-01  4.956405e-01   -6.955053
#> g6  -0.1761358 0.5960219  -1.330530  1.838666e-01  5.253331e-01   -7.115052
#> g20  0.1474773 0.6529736   1.126593  2.603806e-01  5.422614e-01   -7.365100
#> g18  0.1405957 0.5648386   1.057823  2.905765e-01  5.422614e-01   -7.440103
#> g13 -0.1391035 0.5997454  -1.048037  2.950581e-01  5.422614e-01   -7.450394

## plot estimated parameters
plotParameters(pfit)


```
