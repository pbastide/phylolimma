#'
#' @name phyloDE-package
#' @aliases phyloDE
#' @import limma
#' @title The `phyloDE` Package for Differential Expression Analysis
#' @author Paul Bastide, Mélina Gallopin, Arnaud Liehrmann
#' @keywords package
#' @description
#' The `phyloDE` package fits linear model on inter-species gene expression data,
#' combining Phylogenetic Comparative Methods implemented in \code{\link[phylolm]{phylolm}}
#' with moderated statistics tailored for gene expression implemented in \code{\link[limma]{limma}}.
#'
#' With a design matrix that expresses a grouping conditions at the tip of a phylogeny,
#' `phyloDE` can perform Differential Expression analysis.
#'
#' The main function of the package is \code{\link{phylolmFit}},
#' that inherits from the interfaces of both \code{\link[phylolm]{phylolm}}
#' and \code{\link[limma]{limma}}.
#'
"_PACKAGE"
