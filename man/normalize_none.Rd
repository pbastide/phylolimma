% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalize_none}
\alias{normalize_none}
\title{Normalize RNASeq count data}
\usage{
normalize_none(countMatrix, normalisationFactor, dataTransformation)
}
\arguments{
\item{countMatrix}{The RNASeq count matrix. Rows and columns should be named.}

\item{normalisationFactor}{Normalization factors to scale the raw library sizes, as computed e.g. by \code{\link[edgeR]{calcNormFactors}}.}

\item{dataTransformation}{one of "log2", "asin(sqrt)" or "sqrt." See details.}
}
\value{
A matrix of normalized count, with the same dimensions as \code{countMatrix}.
}
\description{
Apply standard CPM, with no length normalization.
}
\keyword{internal}
