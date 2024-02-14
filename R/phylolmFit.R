#' @title Phylogenetic COmparative Method using LIMMA
#'
#' @description
#' This function applies \code{\link[limma]{lmFit}} to the normalized data,
#' in order to take the phylogeny into account.
#' TODO: explain more.
#'
#' @param object 	A matrix data object containing normalized expression values,
#' with rows corresponding to genes and columns to samples (species).
#' @param design the design matrix of the experiment,
#' with rows corresponding to samples and columns to coefficients to be estimated.
#' Defaults to the unit vector (intercept).
#' @param phy an object of class \code{\link[ape]{phylo}},
#' with tips having the same names as the columns of \code{object}.
#' @param model the phylogenetic model used to correct for the phylogeny.
#' Must be one of "BM", "lambda", "OUfixedRoot", "OUrandomRoot" or "delta".
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param measurement_error a logical value indicating whether there is measurement error.
#' Default to \code{TRUE}.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param use_consensus If \code{TRUE}, one consensus tree is used to represent the correlation structure. see \code{\link{phylogeneticCorrelations}}.
#' If \code{FALSE}, each gene will use its own model parameters and will have its own correlation structure accordingly.
#' Default to TRUE.
#' @param consensus_tree If not \code{NULL}, the consensus tree containing the correlation structure,
#' result of function \code{\link{phylogeneticCorrelations}}.
#' If provided, arguments \code{phy}, \code{model} and \code{measurement_error} will be ignored.
#' @param ddf_method the method for the computation of the degrees of freedom of the t statistics (before moderation).
#' Default to \code{ddf_method="Satterthwaite"}.
#' If \code{ddf_method="Species"}, then the number of species is taken for the
#' computation of the degrees of freedom,
#' while if \code{ddf_method="Samples"} the total number of individuals is used.
#' @param REML Use REML (default) or ML for estimating the parameters.
#' @param ... further parameters to be passed
#' to \code{\link[limma]{lmFit}} or \code{\link[phylolm]{phylolm}}.
#'
#' @return An object of class \code{\link[limma]{MArrayLM-class}},
#' with list components \code{coefficients}, \code{stdev.unscaled},
#' \code{sigma} and \code{df.residual}.
#' These quantities take the phylogenetic model into account.
#' The object can be passed to \code{\link[limma]{eBayes}}.
#'
#' @details
#' The default bounds on the phylogenetic parameters are the same as in
#' \code{\link[phylolm]{phylolm}}, except for the \code{alpha} parameter of the OU.
#'
#'
#' @importFrom methods new
#'
#' @export
#'
phylolmFit <- function(object, design = NULL, phy,
                       model = c("BM", "lambda", "OUfixedRoot", "OUrandomRoot", "delta"),
                       measurement_error = FALSE,
                       use_consensus = TRUE,
                       consensus_tree = NULL,
                       ddf_method = c("Satterthwaite", "Species", "Samples"),
                       REML = TRUE, ...) {

  ##################################################################################################
  ## Check unused parameters
  dot_args <- dots(...)
  if ("ndups" %in% names(dot_args) && dot_args$ndups != 1) stop("'ndups' can only be '1' in 'phylolmFit' (for now).")
  if ("spacing" %in% names(dot_args) && dot_args$spacing != 1) stop("'spacing' can only be '1' in 'phylolmFit' (for now).")
  if ("weights" %in% names(dot_args) && !is.null(dot_args$weights)) stop("'weights' can only be 'null' in 'phylolmFit' (for now).")
  if ("method" %in% names(dot_args) && dot_args$method != "ls") stop("'method' can only be 'ls' in 'phylolmFit' (for now).")
  if ("correlation" %in% names(dot_args)) stop("'correlation' is not used in 'phylolmFit' (for now).")
  if ("block" %in% names(dot_args) && !is.null(dot_args$block)) stop("'block' can only be 'null' in 'phylolmFit' (for now).")

  ## Expression Matrix
  if (!is.matrix(object)) stop("'object' must be a matrix.")
  y <- limma::getEAWP(object)
  if (!nrow(y$exprs)) stop("expression matrix has zero rows")

  ##	Check design matrix
  if(is.null(design)) design <- y$design
  if(is.null(design)) {
    design <- matrix(1, ncol(y$exprs), 1)
    rownames(design) <- phy$tip.label
  } else {
    design <- as.matrix(design)
    if(mode(design) != "numeric") stop("design must be a numeric matrix")
    if(nrow(design) != ncol(y$exprs)) stop("row dimension of design doesn't match column dimension of data object")
  }
  ne <- limma::nonEstimable(design)
  if(!is.null(ne)) stop("Coefficients not estimable: ", paste(ne, collapse = " "), "\n")

  ## phylo model
  # if (model != "BM") stop("'modelphy' can only be 'BM' (for now).")
  model <- match.arg(model)

  ## tree
  if (!inherits(phy, "phylo")) stop("object 'phy' must be of class 'phylo'.")
  y_data <- checkParamMatrix(y$exprs, "expression matrix", phy)
  design <- checkParamMatrix(design, "design matrix", phy, transpose = TRUE)

  ## ddf
  ddf_method <- match.arg(ddf_method)

  ##################################################################################################
  ## Consensus tree

  if (!use_consensus && !is.null(consensus_tree)) stop("You set `use_consensus=FALSE`, but provided a consensus tree. Please either set `use_consensus=TRUE` or `consensus_tree=NULL`.")

  if (model == "BM" && !measurement_error) use_consensus <- TRUE

  if (use_consensus) {
    if (!is.null(consensus_tree)) {

      check.consensus_tree(consensus_tree, model, measurement_error)

    } else {
      consensus_tree <- phylogeneticCorrelations(object = object, design = design, phy = phy,
                                                 model = model,
                                                 measurement_error = measurement_error,
                                                 REML = REML,
                                                 ddf_method = ddf_method,
                                                 weights = NULL, ...)
    }

    C_tree_params <- get_chol_tree(y_data, design, consensus_tree$tree, model = "BM", measurement_error = FALSE, REML, ddf_method, ...) ## BM on the consensus tree
    C_tree <- C_tree_params$C_tree
    C_tree_params$optpar <- consensus_tree$alpha
    C_tree_params$lambda_error <- consensus_tree$lambda_error

    ddf_fits <- consensus_tree$ddf

  } else {
    ## one phylo model per gene
    C_tree_params <- get_chol_tree(y_data,  design, phy, model, measurement_error, REML, ddf_method, ...)
    C_tree <- C_tree_params$C_tree

    ddf_fits <- C_tree_params$ddf
  }

  ##################################################################################################

  ## Transform design and data
  design_trans <- transform_design_tree(C_tree, design)
  y_trans <- t(transform_data_tree(C_tree, y_data))

  ## Apply lmFit
  resLmFit <- lmFitLimma(y_trans, design_trans, ...)


  ## Format
  if (use_consensus) {
    resFitFormat <- new("PhyloMArrayLM",
                        list(coefficients = resLmFit$coefficients,
                             sigma = resLmFit$sigma,
                             stdev.unscaled = resLmFit$stdev.unscaled,
                             df.residual = resLmFit$df.residual,
                             Amean = resLmFit$Amean,
                             qr = resLmFit$qr))
  } else {
    resFitFormat <- new("PhyloMArrayLM",
                        list(coefficients = do.call(rbind, resLmFit["coefficients", ]),
                             sigma = do.call(c, resLmFit["sigma", ]),
                             stdev.unscaled = do.call(rbind, resLmFit["stdev.unscaled", ]),
                             df.residual = do.call(c, resLmFit["df.residual", ]),
                             Amean = do.call(c, resLmFit["Amean", ]),
                             qr = resLmFit["qr", ]))
  }

  resFitFormat$df.residual <- ddf_fits

  ## Add slots
  resFitFormat$phy <- phy
  resFitFormat$modelphy <- model
  resFitFormat$measurement_error <- measurement_error
  resFitFormat$phy_trans <- C_tree_params$trans_tree
  resFitFormat$C_tree <- C_tree_params$C_tree
  resFitFormat$optpar <- C_tree_params$optpar
  resFitFormat$lambda_error <- C_tree_params$lambda_error
  resFitFormat$sigma2_phy <- C_tree_params$sigma2_phy
  resFitFormat$sigma2_error <- C_tree_params$sigma2_error
  resFitFormat$REML <- REML
  if (use_consensus) resFitFormat$consensus_tree <- consensus_tree
  resFitFormat$use_consensus <- use_consensus

  ## Result
  return(resFitFormat)
}

#' @title Capture dot arguments
#'
#' @description http://adv-r.had.co.nz/Computing-on-the-language.html#capturing-dots
#'
#' @param ... dots arguments to be captured
#'
#' @return a named list of the arguments in ...
#'
#' @keywords internal
#'
dots <- function(...) {
  eval(substitute(alist(...)))
}

#' @title Fit using limma
#'
#' @description Fit using limma
#'
#'
#' @param y_trans A matrix data object containing normalized
#' and phylogeny transformed expression values,
#' with rows corresponding to genes and columns to samples (species).
#' @param design_trans the phylogeny transformed design matrix of the experiment,
#' with rows corresponding to samples and columns to coefficients to be estimated.
#' Defaults to the unit vector (intercept).
#' @param ... further parameters to be passed to \code{\link[limma]{lmFit}}.
#'
#' @return A list with all the results.
#'
#' @keywords internal
#'
lmFitLimma <- function(y_trans, design_trans, ...) {
  if (!is.list(design_trans)) {
    return(limma::lmFit(object = y_trans, design = design_trans, ...))
  }
  argsfun <- c(ndups = 1, spacing = 1, block = NULL, weights = NULL, method = "ls", dots(...))
  all_res_fit <- mapply(limma::lmFit,
                        object = lapply(seq_len(nrow(y_trans)), function(i) y_trans[i,]),
                        design = design_trans,
                        MoreArgs = argsfun)
  all_res_fit
}

#' @title Get Tree Normalizing Inverse Cholesky
#'
#' @description
#' Compute the whitening cholesky matrix.
#'
#' @param y_data 	A matrix data object containing normalized expression values,
#' with rows corresponding to genes and columns to samples (species).
#' @inheritParams phylolmFit
#'
#' @return The (list of) cholesky matrix of the tree structure.
#'
#' @keywords internal
#'
get_chol_tree <- function(y_data, design, phy, model, measurement_error, REML, ddf_method, ...) {
  if (!measurement_error && model == "BM") { ## Easy case, not fit necessary
    C_tree <- ape::vcv(phy)
    C_tree_chol <- chol(C_tree)
    return(list(C_tree = C_tree_chol,
                trans_tree = phy,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = NA,
                sigma2_error = 0))
  } else {
    C_tree_chol_and_params <- apply(y_data, 1,
                                    get_C_tree, design, phy, model, measurement_error, REML, ddf_method, ...)
    C_tree_chol_and_params <- format_list(C_tree_chol_and_params)
    return(C_tree_chol_and_params)
  }
}

format_list <- function(C_tree_chol_and_params) {
  return(list(C_tree = lapply(C_tree_chol_and_params, function(x) x$C_tree),
              trans_tree = lapply(C_tree_chol_and_params, function(x) x$trans_tree),
              optpar = sapply(C_tree_chol_and_params, function(x) x$optpar),
              lambda_error = sapply(C_tree_chol_and_params, function(x) x$lambda_error),
              sigma2_phy = sapply(C_tree_chol_and_params, function(x) x$sigma2_phy),
              sigma2_error = sapply(C_tree_chol_and_params, function(x) x$sigma2_error),
              ddf = sapply(C_tree_chol_and_params, function(x) x$ddf)))
}

#' @title Get Tree Normalizing Inverse Cholesky
#'
#' @description
#' Compute the whitening cholesky matrix.
#'
#' @param y 	A vector data containing normalized expression values for one gene.
#' @inheritParams phylolmFit
#'
#' @return The cholesky matrix of the tree structure.
#'
#' @keywords internal
#'
get_C_tree <- function(y, design, phy, model, measurement_error, REML, ddf_method, ...) {
  trans_tree_params <- transform_tree_phylolm(y, design, phy, model, measurement_error, REML, ddf_method, ...)
  tree_model <- trans_tree_params$tree_model
  C_tree <- ape::vcv(tree_model)
  C_tree_chol <- chol(C_tree)
  return(list(C_tree = C_tree_chol,
              trans_tree = tree_model,
              optpar = trans_tree_params$optpar,
              lambda_error = trans_tree_params$lambda_error,
              sigma2_phy = trans_tree_params$sigma2_phy,
              sigma2_error = trans_tree_params$sigma2_error,
              ddf = trans_tree_params$ddf))
}

#' @title Get transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_phylolm <- function(y, design, phy, model, measurement_error, REML, ddf_method, ...) {
  if (model == "BM" && !measurement_error) return(phy) # no transformation needed
  data_phylolm <- as.data.frame(cbind(y, design))
  colnames(data_phylolm)[1] <- "expr"
  alpha_bounds <- getBoundsSelectionStrength(phy)
  min_error <- getMinError(phy)
  lower_bounds <- get_lower_bounds(alpha_bounds, min_error, ...)
  upper_bounds <- get_upper_bounds(alpha_bounds, min_error, ...)
  starting_values <- get_starting_values(alpha_bounds, ...)
  dots_args <- get_dots_args(...)
  tmp_fun <- function(...) {
    return(withCallingHandlers(phylolm::phylolm(expr ~ -1 + .,
                                                data = data_phylolm, phy = phy, model = model,
                                                measurement_error = measurement_error,
                                                lower.bound = lower_bounds,
                                                upper.bound = upper_bounds,
                                                starting.value = starting_values,
                                                REML = REML,
                                                ...),
                               warning = function(cond) {
                                 if (grepl(pattern="upper/lower", x = conditionMessage(cond)) && "warning" %in% class(cond)) invokeRestart("muffleWarning")
                               }))
  }
  fplm <- do.call(tmp_fun, dots_args)
  phy_trans_params <- switch(model,
                             BM = transform_tree_model_BM(phy, fplm, measurement_error),
                             lambda = transform_tree_model_lambda(phy, fplm, measurement_error),
                             OUfixedRoot = transform_tree_model_OUfixedRoot(phy, fplm, measurement_error),
                             OUrandomRoot = transform_tree_model_OUrandomRoot(phy, fplm, measurement_error),
                             delta = transform_tree_model_delta(phy, fplm, measurement_error))
  phy_trans_params$ddf <- get_ddf(ddf_method)(fplm, phy)
  return(phy_trans_params)
}

#' @title Get lambda transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_lambda <- function(phy, phyfit, measurement_error) {
  if (measurement_error) stop("Measurement error is not allowed with lambda model.")
  return(list(tree_model = phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = phyfit$optpar))$tree,
              optpar = NA,
              lambda_error = 1,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = 0))
}

#' @title Get BM transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_BM <- function(phy, phyfit, measurement_error) {
  if (!measurement_error) stop("Measurement error should be true here.")
  lambda_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tree_height(phy))
  return(list(tree_model = phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_error))$tree,
              optpar = NA,
              lambda_error = lambda_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Get OU transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_OUfixedRoot <- function(phy, phyfit, measurement_error) {
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = phyfit$optpar))$tree
  if (!measurement_error) {
    return(list(tree_model = tree_model,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = phyfit$sigma2,
                sigma2_error = 0))
  }
  tilde_t <- tree_height(tree_model) / (2 * phyfit$optpar)
  lambda_ou_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_ou_error))$tree
  return(list(tree_model = tree_model,
              optpar = phyfit$optpar,
              lambda_error = lambda_ou_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Get OU transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_OUrandomRoot <- function(phy, phyfit, measurement_error) {
  tree_model <- phylolm::transf.branch.lengths(phy, "OUrandomRoot", parameters = list(alpha = phyfit$optpar))$tree
  tree_model$root.edge <- 0
  if (!measurement_error) {
    return(list(tree_model = tree_model,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = phyfit$sigma2,
                sigma2_error = 0))
  }
  tilde_t <- tree_height(tree_model) / (2 * phyfit$optpar)
  lambda_ou_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_ou_error))$tree
  return(list(tree_model = tree_model,
              optpar = phyfit$optpar,
              lambda_error = lambda_ou_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Get delta transformed tree
#'
#' @description
#' Compute the transformed tree using \code{\link[phylolm]{transf.branch.lengths}}.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
transform_tree_model_delta <- function(phy, phyfit, measurement_error) {
  tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = phyfit$optpar))$tree
  if (!measurement_error) {
    return(list(tree_model = tree_model,
                optpar = NA,
                lambda_error = 1,
                sigma2_phy = phyfit$sigma2,
                sigma2_error = 0))
  }
  tilde_t <- tree_height(tree_model)
  lambda_delta_error <- get_lambda_error(phyfit$sigma2, phyfit$sigma2_error, tilde_t)
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_delta_error))$tree
  return(list(tree_model = tree_model,
              optpar = phyfit$optpar,
              lambda_error = lambda_delta_error,
              sigma2_phy = phyfit$sigma2,
              sigma2_error = phyfit$sigma2_error))
}

#' @title Transform design matrix
#'
#' @description
#' Multiply by inverse cholesky to whiten the data
#'
#' @param C_tree Cholesky of the tree structure obtained through \code{\link{get_chol_tree}}
#' @param design design matrix
#'
#' @return transformed design matrix
#'
#' @keywords internal
#'
transform_design_tree <- function(C_tree, design) {
  if (!is.list(C_tree)) return(transform_design_one_tree(C_tree, design))
  return(lapply(C_tree, transform_design_one_tree, design))
}

transform_design_one_tree <- function(C_tree, design, transpose = FALSE) {
  res <- backsolve(C_tree, design, transpose = TRUE)
  colnames(res) <- colnames(design)
  rownames(res) <- rownames(design)
  return(res)
}

#' @title Transform data matrix
#'
#' @description
#' Multiply by inverse cholesky to whiten the data
#'
#' @param C_tree Cholesky of the tree structure obtained through \code{\link{get_chol_tree}}
#' @param y_data data matrix
#'
#' @return transformed data matrix
#'
#' @keywords internal
#'
transform_data_tree <- function(C_tree, y_data) {
  if (!is.list(C_tree)) return(transform_design_one_tree(C_tree, t(y_data)))
  return(mapply(transform_design_one_tree,
                C_tree,
                lapply(seq_len(nrow(y_data)), function(i) y_data[i,])))
}

setGeneric("log_likelihood", function(object) standardGeneric("log_likelihood"))
setMethod("log_likelihood", "MArrayLM", function(object) NULL)
setMethod("log_likelihood", "PhyloMArrayLM", function(object) log_likelihood_internal(object))

log_likelihood_internal <- function (object) {
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
