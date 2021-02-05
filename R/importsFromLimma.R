#' @title Class PhyloMArrayLM
#'
#' @description
#' Extension of class \code{\link[limma]{MArrayLM-class}} from package \code{limma}
#' for phylogenetically correlated regressions of gene-wise linear models.
#'
#' It is a simple list-based S4 class.
#' Objects are normally created by \code{\link{phylolmFit}}.
#' Additional components are added by \code{\link[limma]{eBayes}}.
#'
#' @section Components:
#' \code{PhyloMArrayLM} objects do not contain any slots (apart from .Data)
#' but they should contain the same list components than \code{\link[limma]{MArrayLM-class}}.
#'
#' In addition, they contain the following tree-specific components:
#'
#' \describe{
#' \item{\code{phy}:}{The phylogenetic tree used for the regression, of class \code{\link[ape]{phylo}}.}
#' \item{\code{modelphy}:}{The phylogenetic model of trait evolution, argument call in \code{\link{phylolmFit}}.}
#' \item{\code{measurement_error}:}{Boolean, TRUE if there is additional measurement error, argument call in \code{\link{phylolmFit}}.}
#' \item{\code{phy_trans}:}{List of the transformed phylogenetic trees obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' \item{\code{optpar}:}{Vector of the \code{optpar} obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' \item{\code{lambda_error}:}{Vector of the \code{optpar} parameters obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' \item{\code{sigma2_phy}:}{Vector of the phylogenetic \code{sigma2} parameters obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' \item{\code{sigma2_error}:}{Vector of the \code{sigma2_error} parameters obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' }
#'
#' @seealso \code{\link[limma]{MArrayLM-class}}, \code{\link{phylolmFit}}.
#'
#' @export
#' @import methods
#' @importClassesFrom limma MArrayLM
.PhyloMArrayLM <- setClass("PhyloMArrayLM", representation("list"), contains="MArrayLM")


#' @title Empirical Bayes Statistics for Differential Expression
#'
#' @description
#' This function applies \code{\link[limma]{eBayes}} to the result of function
#' \code{\link{phylolmFit}}.
#'
#' @param fit a \code{\linkS4class{PhyloMArrayLM}} object, fitted using \code{\link[limma]{eBayes}}.
#' @param ... further parameters to be passed to \code{\link[limma]{eBayes}}.
#'
#' @export
#'
eBayes <- function(fit, ...) {
  return(limma::eBayes(fit, ...))
}
