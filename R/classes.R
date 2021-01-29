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
#' \item{\code{TODO}:}{TODO}
#' }
#'
#' @seealso \code{\link[limma]{MArrayLM-class}}, \code{\link{phylolmFit}}.
#'
#' @export
#' @import methods
#' @importClassesFrom limma MArrayLM
.PhyloMArrayLM <- setClass("PhyloMArrayLM", representation("list"), contains="MArrayLM")
