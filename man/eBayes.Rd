% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importsFromLimma.R
\name{eBayes}
\alias{eBayes}
\alias{treat}
\title{Empirical Bayes Statistics for Differential Expression}
\usage{
eBayes(fit, ...)

treat(fit, ...)
}
\arguments{
\item{fit}{a \code{\linkS4class{PhyloMArrayLM}} object, fitted using \code{\link{phylolmFit}}.}

\item{...}{further parameters to be passed to \code{\link[limma]{eBayes}} or \code{\link[limma]{treat}}.}
}
\description{
Apply \code{\link[limma]{eBayes}} to the result of function
\code{\link{phylolmFit}}.

Function \code{\link[limma]{treat}} is not supported yet for a \code{\linkS4class{PhyloMArrayLM}},
and will throw an error.
}
