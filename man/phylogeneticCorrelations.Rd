% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylogeneticCorrelations.R
\name{phylogeneticCorrelations}
\alias{phylogeneticCorrelations}
\title{Phylogenetic COmparative Method using LIMMA}
\usage{
phylogeneticCorrelations(
  object,
  design = NULL,
  phy,
  col_species = NULL,
  model = c("BM", "lambda", "OUfixedRoot", "OUrandomRoot", "delta"),
  measurement_error = TRUE,
  trim = 0.15,
  weights = NULL,
  REML = TRUE,
  ddf_method = c("Satterthwaite", "Species", "Samples"),
  ncores = 1,
  ...
)
}
\arguments{
\item{object}{A matrix data object containing normalized expression values,
with rows corresponding to genes and columns to samples (species).}

\item{design}{the design matrix of the experiment,
with rows corresponding to samples and columns to coefficients to be estimated.
Defaults to the unit vector (intercept).}

\item{phy}{an object of class \code{\link[ape]{phylo}}.
It must be either a tree with tips having the same names as the columns of \code{object} (including replicates),
or a tree such that tip labels match with species names in `col_species`.}

\item{col_species}{a character vector with same length as columns in the expression matrix,
specifying the species for the corresponding column. If left `NULL`, an automatic parsing of species names with sample ids is attempted.}

\item{model}{the phylogenetic model used to correct for the phylogeny.
Must be one of "BM", "lambda", "OUfixedRoot", "OUrandomRoot" or "delta".
See \code{\link[phylolm]{phylolm}} for more details.}

\item{measurement_error}{a logical value indicating whether there is measurement error.
Default to \code{TRUE}.
See \code{\link[phylolm]{phylolm}} for more details.}

\item{trim}{the fraction of observations to be trimmed from each end when computing the trimmed mean. Default to 0.15, as in \code{\link[limma]{duplicateCorrelation}}.}

\item{weights}{a named vector or matrix with weights to be applied on the measurement error term.
See \code{\link[phylolm]{phylolm}} for more details.}

\item{REML}{Use REML (default) or ML for estimating the parameters.}

\item{ddf_method}{the method for the computation of the degrees of freedom of the t statistics (before moderation).
Default to \code{ddf_method="Satterthwaite"}.
If \code{ddf_method="Species"}, then the number of species is taken for the
computation of the degrees of freedom,
while if \code{ddf_method="Samples"} the total number of individuals is used.}

\item{ncores}{number of cores to use for parallel computation. Default to 1 (no parallel computation).}

\item{...}{further parameters to be passed
to \code{\link[limma]{lmFit}}.}
}
\value{
An object of class \code{TransTree-class},
with list components:
\itemize{
\item \code{tree} the transformed consensus tree
\item \code{params} the associated consensus parameters.
}
This consensus tree defines a correlation structure, and can be passed on to \code{\link[limma]{lmFit}}.
}
\description{
This function applies \code{\link[limma]{lmFit}} to the normalized data,
in order to take the phylogeny into account.
TODO: explain more.
}
\seealso{
\code{\link[limma]{lmFit}}, \code{\link[phylolm]{phylolm}}, \code{\link[limma]{duplicateCorrelation}}
}
