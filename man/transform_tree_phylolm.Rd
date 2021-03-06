% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylolmFit.R
\name{transform_tree_phylolm}
\alias{transform_tree_phylolm}
\title{Get transformed tree}
\usage{
transform_tree_phylolm(y, design, phy, model, measurement_error, ...)
}
\arguments{
\item{y}{A vector data containing normalized expression values for one gene.}

\item{design}{the design matrix of the experiment,
with rows corresponding to samples and columns to coefficients to be estimated.
Defaults to the unit vector (intercept).}

\item{phy}{an object of class \code{\link[ape]{phylo}},
with tips having the same names as the columns of \code{object}.}

\item{model}{the phylogenetic model used to correct for the phylogeny.
Must be one of "BM", "lambda", "OUfixedRoot" or "delta".
See \code{\link[phylolm]{phylolm}} for more details.}

\item{measurement_error}{a logical value indicating whether there is measurement error.
See \code{\link[phylolm]{phylolm}} for more details.}

\item{...}{further parameters to be passed
to \code{\link[limma]{lmFit}} or \code{\link[phylolm]{phylolm}}.}
}
\value{
The transformed tree.
}
\description{
Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
}
\keyword{internal}
