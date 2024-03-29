% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylogeneticCorrelations.R
\name{get_consensus_tree_OUrandomRoot}
\alias{get_consensus_tree_OUrandomRoot}
\title{Get OU transformed tree}
\usage{
get_consensus_tree_OUrandomRoot(
  phy,
  all_phyfit,
  measurement_error,
  trim,
  median = FALSE
)
}
\arguments{
\item{phy}{an object of class \code{\link[ape]{phylo}}.
It must be either a tree with tips having the same names as the columns of \code{object} (including replicates),
or a tree such that tip labels match with species names in `col_species`.}

\item{measurement_error}{a logical value indicating whether there is measurement error.
Default to \code{TRUE}.
See \code{\link[phylolm]{phylolm}} for more details.}

\item{median}{if TRUE, the \code{pracma::geo_median} function is used to take
the geometric median.}
}
\value{
The transformed tree.
}
\description{
Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
}
\keyword{internal}
