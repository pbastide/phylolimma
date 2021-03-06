% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalize_RPKM}
\alias{normalize_RPKM}
\title{Normalize RNASeq count data}
\usage{
normalize_RPKM(
  countMatrix,
  lengthMatrix,
  normalisationFactor,
  dataTransformation
)
}
\arguments{
\item{countMatrix}{The RNASeq count matrix. Rows and columns should be named.}

\item{lengthMatrix}{The associated length matrix. Should have the same dimensions as \code{countMatrix}, with the same names.}

\item{normalisationFactor}{Normalization factors to scale the raw library sizes, as computed e.g. by \code{\link[edgeR]{calcNormFactors}}.}

\item{dataTransformation}{one of "log2", "asin(sqrt)" or "sqrt." See details.}
}
\value{
A matrix of normalized count, with the same dimensions as \code{countMatrix}.
}
\description{
Apply RPKM length normalization.
}
\keyword{internal}
