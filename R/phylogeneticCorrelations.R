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
#' Must be one of "BM", "lambda", "OUfixedRoot" or "delta".
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param measurement_error a logical value indicating whether there is measurement error.
#' Default to \code{TRUE}.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param trim the fraction of observations to be trimmed from each end when computing the trimmed mean. Default to 0.15, as in \code{\link[limma]{duplicateCorrelation}}.
#' @param weights a named vector or matrix with weights to be applied on the measurement error term.
#' See \code{\link[phylolm]{phylolm}} for more details.
#' @param ... further parameters to be passed
#' to \code{\link[limma]{lmFit}} or \code{\link[phylolm]{phylolm}}.
#'
#' @return An object of class \code{TransTree-class},
#' with list components:
#' \itemize{
#' \item \code{tree} the transformed consensus tree
#' \item \code{params} the associated consensus parameters.
#' }
#' This consensus tree define a correlation structure, and can be passed on to \code{\link[limma]{lmFit}}.
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[phylolm]{phylolm}}, \code{\link[limma]{duplicateCorrelation}}
#'
#' @importFrom methods new
#'
#' @export
#'
#' @importFrom graphics lines title
#' @importFrom stats approxfun lowess model.matrix
#'
phylogeneticCorrelations <- function(object, design = NULL, phy,
                                     model = c("BM", "lambda", "OUfixedRoot", "delta"),
                                     measurement_error = TRUE,
                                     trim = 0.15, weights = NULL, ...) {

  ##################################################################################################
  ## Checks

  ## Expression Matrix
  if (!is.matrix(object)) stop("'object' must be a matrix.")
  y <- limma::getEAWP(object)
  if (!nrow(y$exprs)) stop("expression matrix has zero rows")

  #	Check weights
  if(!is.null(weights)) {
    stop("weights are not allowed with the phylogenetic regression (yet).")
    # message("'weights' will be used in the independent errors.")
    # weights <- limma::asMatrixWeights(weights, dim(y))
    # weights[weights <= 0] <- NA
    # y[!is.finite(weights)] <- NA
  }

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
  if(!is.null(ne)) cat("Coefficients not estimable:", paste(ne, collapse = " "), "\n")

  ## phylo model
  model <- match.arg(model)

  ## tree
  if (!inherits(phy, "phylo")) stop("object 'phy' must be of class 'phylo'.")
  y_data <- checkParamMatrix(y$exprs, "expression matrix", phy)
  design <- checkParamMatrix(design, "design matrix", phy, transpose = TRUE)

  ##################################################################################################

  tree_model <- get_consensus_tree(y_data, design, phy, model, measurement_error, weights, trim, ...)

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
get_consensus_tree <- function(y_data, design, phy, model, measurement_error, weights, trim, ...) {
  if(!is.null(weights)) stop("weights are not allowed with the phylogenetic regression (yet).")

  if (model == "BM" && !measurement_error) # no parameter to estimate
    return(list(tree = phy,
                params = list(model = "BM",
                              measurement_error = FALSE)))

  all_fits <- list()
  for (i in 1:nrow(y_data)) {
    y <- y_data[i, ]
    # w <- weights[i, ]

    data_phylolm <- as.data.frame(cbind(y, design))
    colnames(data_phylolm)[1] <- "expr"
    lower_bounds <- get_lower_bounds(...)
    dots_args <- get_dots_args(...)

    tmp_fun <- function(...) {
      return(phylolm::phylolm(expr ~ -1 + ., data = data_phylolm, phy = phy, model = model,
                              measurement_error = measurement_error, lower.bound = lower_bounds,
                              ...))
                              # error_weight = weights, ...))
    }
    all_fits[[i]] <- do.call(tmp_fun, dots_args)

  }

  params <- switch(model,
                   BM = get_consensus_tree_BM(phy, all_fits, measurement_error, trim),
                   lambda = get_consensus_tree_lambda(phy, all_fits, measurement_error, trim),
                   OUfixedRoot = get_consensus_tree_OUfixedRoot(phy, all_fits, measurement_error, trim),
                   delta = get_consensus_tree_delta(phy, all_fits, measurement_error, trim))
  return(params)
}

#' @title Get lambda transformed tree
#'
#' @description
#' Compute the consensus parameters.
#'
#' @inheritParams get_C_tree
#'
#' @return The transformed tree.
#'
#' @keywords internal
#'
get_consensus_tree_lambda <- function(phy, all_phyfit, measurement_error, trim) {
  if (measurement_error) stop("Measurement error is not allowed with lambda model.")

  all_lambdas <- sapply(all_phyfit, function(x) x$optpar)
  all_lambdas_transform <- atanh(pmax(-1, all_lambdas))
  lambda_mean <- tanh(mean(all_lambdas_transform, trim = trim, na.rm = TRUE))

  return(list(tree = phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_mean))$tree,
              params = list(model = "lambda",
                            measurement_error = measurement_error,
                            lambda = lambda_mean)))
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
get_consensus_tree_BM <- function(phy, all_phyfit, measurement_error, trim) {

  if (!measurement_error) stop("Measurement error should be true here.")

  all_lambda_error <- sapply(all_phyfit, function(phyfit) phyfit$sigma2 / (phyfit$sigma2_error / max(ape::vcv(phy)) + phyfit$sigma2))

  all_lambdas_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_mean <- tanh(mean(all_lambdas_transform, trim = trim, na.rm = TRUE))

  return(list(tree = phylolm::transf.branch.lengths(phy, "lambda", parameters = list(lambda = lambda_mean))$tree,
              params = list(model = "BM",
                            measurement_error = measurement_error,
                            lambda_error = lambda_mean)))
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
get_consensus_tree_OUfixedRoot <- function(phy, all_phyfit, measurement_error, trim) {

  all_alphas <- sapply(all_phyfit, function(x) x$optpar)
  all_alphas_transform <- log(all_alphas)
  alpha_mean <- exp(mean(all_alphas_transform, trim = trim, na.rm = TRUE))

  if (!measurement_error) {

    return(list(tree = phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_mean))$tree,
                params = list(model = "OUfixedRoot",
                              measurement_error = measurement_error,
                              alpha = alpha_mean)))
  }

  get_lambda_error <- function(phyfit) {
    tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
                                                 parameters = list(alpha = phyfit$optpar))$tree
    tilde_t <- max(ape::vcv(tree_model)) / (2 * phyfit$optpar)
    lambda_ou_error <- phyfit$sigma2 * tilde_t / (phyfit$sigma2_error + phyfit$sigma2 * tilde_t)
    return(lambda_ou_error)
  }

  ## consensus lambda error
  all_lambda_error <- sapply(all_phyfit, get_lambda_error)
  all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))

  ## transform tree
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot", parameters = list(alpha = alpha_mean))$tree
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree

  return(list(tree = tree_model,
              params = list(model = "OUfixedRoot",
                            measurement_error = measurement_error,
                            alpha = alpha_mean,
                            lambda_error = lambda_error_mean)))
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
get_consensus_tree_delta <- function(phy, all_phyfit, measurement_error, trim) {

  all_deltas <- sapply(all_phyfit, function(x) x$optpar)
  all_deltas_transform <- log(all_deltas)
  delta_mean <- exp(mean(all_deltas_transform, trim = trim, na.rm = TRUE))

  if (!measurement_error) {

    return(list(tree = phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = delta_mean))$tree,
                params = list(model = "delta",
                              measurement_error = measurement_error,
                              delta = delta_mean)))
  }

  get_lambda_error <- function(phyfit) {
    tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = phyfit$optpar))$tree
    tilde_t <- max(ape::vcv(tree_model))
    lambda_delta_error <- phyfit$sigma2 * tilde_t / (phyfit$sigma2_error + phyfit$sigma2 * tilde_t)
    return(lambda_delta_error)
  }

  all_lambda_error <- sapply(all_phyfit, get_lambda_error)
  all_lambda_error_transform <- atanh(pmax(-1, all_lambda_error))
  lambda_error_mean <- tanh(mean(all_lambda_error_transform, trim = trim, na.rm = TRUE))

  ## transform tree
  tree_model <- phylolm::transf.branch.lengths(phy, "delta", parameters = list(delta = delta_mean))$tree
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_error_mean))$tree

  return(list(tree = tree_model,
              params = list(model = "delta",
                            measurement_error = measurement_error,
                            delta = delta_mean,
                            lambda_error = lambda_error_mean)))
}

check.consensus_tree <- function(consensus_tree, model, measurement_error) {
  if (consensus_tree$params$model != model) {
    stop(paste0("The consensus tree was computed with the ", consensus_tree$params$model, " model, but the model used is ", model, ". Please check which model is the correct one."))
  }
  if (consensus_tree$params$measurement_error != measurement_error) {
    stop(paste0("The consensus tree was computed with `measurement_error=", consensus_tree$params$measurement_error, "`, but you specified `measurement_error=", measurement_error, "`. Please check which model is the correct one."))
  }
}
