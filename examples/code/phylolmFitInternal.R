# #' @title Internal phylolm fit
# #'
# #' @description
# #' Lightweight lm fit assuming that the data is already formatted.
# #'
# #' @param X the regression matrix. Should be consistent with the tree and data.
# #' @param y the response vector. Should be in the same order as the tip labels of the tree.
# #' @param phy the phylogenetic tree. Must be ultrametric, and in preorder.
# #' @param method one of "ML" or "REML" (the default).
# #' @param dof_method one of "Satterthwaite" (default), "Species" or "Samples"
# #' @inheritParams phylolm::phylolm
# #'
# #' @return fit
# #'
# #' @export
# #'
#
# phylolmFitInternal <- function(X, y, phy,
#                                model = c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"),
#                                lower.bound = NULL, upper.bound = NULL, starting.value = NULL,
#                                measurement_error = FALSE,
#                                method = c("REML", "ML")) {
#
#
#   n <- nrow(X)
#   d <- ncol(X)
#   reml <- method == "reml"
#   log2pi <- log(2 * pi)
#
#   ## Likelihood
#   if (!reml) {
#     minusLogLik <- function(parameters, X, y, phy, model) {
#       n <- nrow(X)
#       d <- ncol(X)
#       phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
#       comp <- three.point.compute(phytrans, P = y - yhat, Q = y)
#       n2llh <- as.numeric( n * log2pi + n * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2) # -2 log-likelihood
#       return(n2llh / 2)
#     }
#   }
#
#
#
# }
