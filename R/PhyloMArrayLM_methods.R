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
#' \item{\code{lambda_error}:}{Vector of the \code{lambda_error} parameters obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' \item{\code{sigma2_phy}:}{Vector of the phylogenetic \code{sigma2} parameters obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' \item{\code{sigma2_error}:}{Vector of the \code{sigma2_error} parameters obtained from a fit using \code{\link[phylolm]{phylolm}} on each gene.}
#' }
#'
#' @seealso \code{\link[limma]{MArrayLM-class}}, \code{\link{phylolmFit}}.
#'
#' @export
#' @import methods
#' @importClassesFrom limma MArrayLM
#'
setClass("PhyloMArrayLM",
         representation("list"),
         contains = "MArrayLM",
)

#' @title Class ConsensusTreeModel
#'
#' @description
#' A simple list-based S4 class, that contains a tree and associated parameters,
#' obtained through function \code{\link{phylogeneticCorrelations}}.
#'
#' @section Components:
#' \code{ConsensusTreeModel} objects do not contain any slots (apart from .Data)
#' but they should contain the following list components:
#'
#'#' \itemize{
#' \item \code{tree} the transformed consensus tree
#' \item \code{params} the associated consensus parameters.
#' }
#'
#' \describe{
#' \item{\code{tree}:}{the transformed consensus tree used to define the correlation structure.}
#' \item{\code{params}:}{the parameters associated with the tree, including regularized values of the parameters.}
#' }
#'
#' @seealso \code{\link{phylogeneticCorrelations}}, \code{\link{phylolmFit}}, \code{\link{PhyloMArrayLM-class}}.
#'
#' @export
#'
setClass("ConsensusTreeModel",
         representation("list")
)

#' Methods for classes PhyloMArrayLM and ConsensusTreeModel
#'
#' @param object an object of class \code{\link{PhyloMArrayLM-class}} or \code{\link{ConsensusTreeModel-class}}
#'
#' @param consensus if \code{TRUE} (the default), \code{getParameters}
#' returns the consensus parameters as a named vector.
#' Otherwise, it returns the individual fits for all genes.
#'
#' @param ... further parameters to be based to \code{\link{hist}}, including \code{breaks}.
#'
#' @return Depending on the method, a named vector, a data frame, or a plot.
#'
#' @seealso \code{\link{PhyloMArrayLM-class}}, \code{\link{phylolmFit}}
#'
#' @export
#' @docType methods
#' @rdname PhyloMArrayLMMethods
#'
setGeneric("getParameters", function(object, consensus = TRUE) standardGeneric("getParameters"))

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("getParameters", "PhyloMArrayLM", function(object, consensus = TRUE) {
  if (object$modelphy == "BM" && !object$measurement_error) return(NULL)
  if (object$use_consensus) {
    return(getParameters(object$consensus_tree, consensus = consensus))
  } else {
    if (consensus) stop("Fit did not use a consensus tree. Set `consensus = FALSE` to get the individual gene-specific parameters, or re-run fit with consensus tree.")
    params <- data.frame(lambda = object$lambda_error)
    if (object$modelphy == "OUfixedRoot") {
      params$rho <- rhoFromAlpha(object$optpar, max(ape::node.depth.edgelength(object$phy)))
    }
    return(params)
  }
})

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("getParameters", "ConsensusTreeModel", function(object, consensus = TRUE) {
  objpar <- object$param
  if (objpar$model == "BM" && !objpar$measurement_error) return(NULL)
  if (consensus) {
    params <- params_names <- NULL
    if (!is.null(objpar$lambda_error)) {
      params <- objpar$lambda_error
      params_names <- "lambda"
    }
    if (objpar$model == "OUfixedRoot") {
      params <- c(params, rhoFromAlpha(objpar$alpha, objpar$t_original_tree))
      params_names <- c(params_names, "rho")
    }
    names(params) <- params_names
  } else {
    params <- data.frame(lambda = tanh(objpar$atanh_lambda_error))
    if (objpar$model == "OUfixedRoot") {
      params$rho <- tanh(objpar$trans_alpha)
    }
  }
  return(params)
})

#'@importFrom graphics close.screen screen split.screen title hist abline
NULL

#' @rdname PhyloMArrayLMMethods
setGeneric("plotParameters", function(object, ...) standardGeneric("plotParameters"))

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("plotParameters", "PhyloMArrayLM", function(object, ...) {
  params <- NULL
  if (object$use_consensus) params <- getParameters(object, consensus = TRUE)
  plot_params_internal(params, getParameters(object, consensus = FALSE), ...)
})

#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("plotParameters", "ConsensusTreeModel", function(object, ...) {
  plot_params_internal(getParameters(object, consensus = TRUE),
                       getParameters(object, consensus = FALSE), ...)
})

plot_params_internal <- function(params, all_params, ...) {
  ncols <- ncol(all_params)
  scr <- split.screen(c(1, ncols))
  on.exit(close.screen(all.screens = TRUE))
  for (i in seq_len(ncol(all_params))) {
    screen(scr[i])
    hist(all_params[, i], xlim = c(0, 1), xlab = colnames(all_params)[i],
         main = "", ...)
    if (!is.null(params)) abline(v = params[i], lty = "dashed", lwd = 2)
  }
}

#' @rdname PhyloMArrayLMMethods
setMethod("show", "PhyloMArrayLM", function(object) {
  cat(is(object)[[1]], "\n",
      "  Fit on: ", nrow(object$coef), " genes.\n",
      "  Model:  ", object$modelphy, ifelse(object$measurement_error, ", with", ", without"), " measurement error.\n",
      ifelse(object$use_consensus, "  Using", "  Not using"), " a consensus tree.\n",
      sep = ""
  )
})

#' @rdname PhyloMArrayLMMethods
setMethod("show", "ConsensusTreeModel", function(object) {
  cat(is(object)[[1]], "\n",
      "  Consensus tree on: ", length(object$params$atanh_lambda_error), " genes.\n",
      "  Model: ", object$param$model, ifelse(object$params$measurement_error, ", with", ", without"), " measurement error.\n",
      "  Consensus parameters: ",
      sep = ""
  )
  getParametersString(object)
})

getParametersString <- function(object) {
  pp <- getParameters(object)
  cat(names(pp)[1], "=", pp[1])
  if (length(pp) > 1) {
    for (i in 2:length(pp)) {
      cat(",", names(pp)[i], "=", pp[i])
    }
  }
}

#' @rdname PhyloMArrayLMMethods
setGeneric("logLikelihood", function(object) standardGeneric("logLikelihood"))
#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("logLikelihood", "PhyloMArrayLM", function(object) logLikelihood_internal(object))

logLikelihood_internal <- function (object) {
  REML <- object$REML
  sigma_hat <- object$sigma^2
  N <- length(object$phy$tip.label)
  p <-  N - object$df.residual
  sum_res <- sigma_hat * (N - p)
  if (!is.null(object$weights)) stop("A PhyloMArrayLM object cannot have weights.")
  if (!is.list(object$C_tree)) {
    tree_det <- sum(log(diag(object$C_tree)))
  } else {
    tree_det <- sapply(object$C_tree, function(CC) sum(log(diag(CC))))
  }
  N0 <- N
  if (REML) N <- N - p
  val <- 0.5 * (- N * (log(2 * pi) + 1 + log(sum_res) - log(N))) - tree_det
  if (REML) {
    if (object$use_consensus) {
      val <- val - sapply(p, function(pp) sum(log(abs(diag(object$qr$qr)[1L:pp]))))
    } else {
      val <- val - sapply(1:length(p), function(pp) sum(log(abs(diag(object$qr[[pp]]$qr)[1L:p[pp]]))))
    }
  }
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}


#' @rdname PhyloMArrayLMMethods
#' @export
setGeneric("consensusTree", function(object) standardGeneric("consensusTree"))
#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("consensusTree", "PhyloMArrayLM", function(object) {
  if (!object$use_consensus) {
    warning("The fitted object did not use a consensus tree.")
    return(NULL)
  }
  return(consensusTree(object$consensus_tree))
})
#' @rdname PhyloMArrayLMMethods
#' @export
setMethod("consensusTree", "ConsensusTreeModel", function(object) {
  return(object$tree)
})

consensus_tree_internal <- function (object) {
  if (!object$use_consensus) {
    warning("The fitted object did not use a consensus tree.")
    return(NULL)
  }
  return(object$consensus_tree$tree)
}


