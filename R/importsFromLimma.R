#' @title Not Implemented Functions
#'
#' @description
#' Functions
#' \code{\link[limma]{treat}},
#' \code{\link[limma]{decideTests}},
#' \code{\link[limma]{classifyTestsF}}
#' are not supported yet for a \code{\linkS4class{PhyloMArrayLM}},
#' and will throw an error.
#'
#' @param fit a \code{\linkS4class{PhyloMArrayLM}} object, fitted using \code{\link{phylolmFit}}.
#' @param ... further parameters.
#'
#' @export
treat <- function(fit, ...) {
  if (is(fit, "PhyloMArrayLM")) stop("Function `treat` is not supported for an object of class `PhyloMArrayLM`.")
  return(limma::treat(fit, ...))
}

#' @rdname treat
#' @export
decideTests <- function(fit, ...) {
  if (is(fit, "PhyloMArrayLM")) stop("Function `decideTests` is not supported for an object of class `PhyloMArrayLM`.")
  return(limma::decideTests(fit, ...))
}

#' @rdname treat
#' @export
#'
classifyTestsF <- function(fit, ...) {
  if (is(fit, "PhyloMArrayLM")) stop("Function `classifyTestsF` is not supported for an object of class `PhyloMArrayLM`.")
  return(limma::classifyTestsF(fit, ...))
}
