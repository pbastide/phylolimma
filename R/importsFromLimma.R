#' @title Empirical Bayes Statistics for Differential Expression
#'
#' @description
#' Apply \code{\link[limma]{eBayes}} to the result of function
#' \code{\link{phylolmFit}}.
#'
#' Function \code{\link[limma]{treat}} is not supported yet for a \code{\linkS4class{PhyloMArrayLM}},
#' and will throw an error.
#'
#' @param fit a \code{\linkS4class{PhyloMArrayLM}} object, fitted using \code{\link{phylolmFit}}.
#' @param ... further parameters to be passed to \code{\link[limma]{eBayes}} or \code{\link[limma]{treat}}.
#'
#' @export
#'
eBayes <- function(fit, ...) {
  return(limma::eBayes(fit, ...))
}

#'
#' @rdname eBayes
#' @export
#'
treat <- function(fit, ...) {
  if (is(fit, "PhyloMArrayLM")) stop("Function `treat` is not supported for an object of class `PhyloMArrayLM`.")
  return(limma::treat(fit, ...))
}

#' @title Multiple Testing Across Genes and Contrasts
#'
#' @description
#' Function \code{\link[limma]{decideTests}} is not supported yet for a \code{\linkS4class{PhyloMArrayLM}},
#' and will throw an error.
#'
#' @param fit a \code{\linkS4class{PhyloMArrayLM}} object, fitted using \code{\link{phylolmFit}}.
#' @param ... further parameters to be passed to \code{\link[limma]{decideTests}}.
#'
#' @export
#'
decideTests <- function(fit, ...) {
  if (is(fit, "PhyloMArrayLM")) stop("Function `decideTests` is not supported for an object of class `PhyloMArrayLM`.")
  return(limma::decideTests(fit, ...))
}

#' @title Multiple Testing Across Genes and Contrasts
#'
#' @description
#' Function \code{\link[limma]{classifyTestsF}} is not supported yet for a \code{\linkS4class{PhyloMArrayLM}},
#' and will throw an error.
#'
#' @param fit a \code{\linkS4class{PhyloMArrayLM}} object, fitted using \code{\link{phylolmFit}}.
#' @param ... further parameters to be passed to \code{\link[limma]{classifyTestsF}}.
#'
#' @rdname decideTests
#' @export
#'
classifyTestsF <- function(fit, ...) {
  if (is(fit, "PhyloMArrayLM")) stop("Function `classifyTestsF` is not supported for an object of class `PhyloMArrayLM`.")
  return(limma::classifyTestsF(fit, ...))
}
