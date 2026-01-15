# #' @title Initialization of the gamma parameter.
# #'
# #' @description
# #' Initialize the stationary variance of the process.
# #'
# #' @param X the design matrix
# #' @param y the trait vector
# #' @param tree the phylogenetic tree
# #'
# #' @details
# #' Initialize the stationary variance by taking the mad over the centered and
# #' normalized data. This assumes that the process as the tips is at the
# #' stationary state, so that all tips have variance gamma^2 = sigma^2 / (2 alpha).
# #' The data is centered using a (non phylogenetic) robust regression of y
# #' against X, and extracting the fitted values, that are taken to approximate the center.
# #'
# #' @return An initial gamma^2 parameter guess
# #'
# #' @keywords internal
# #'
# initGammaParameter <- function(X, y, tree) {
#   center <- robustbase::lmrob.S(X, y, control = robustbase::lmrob.control())$fitted.values
#   norm_trait <- (y - center)
#   return(mad(norm_trait, center = 0)^2)
# }
#
# #' @title Initialization of the alpha parameter.
# #'
# #' @description
# #' Initialize the stationary variance of the process.
# #'
# #' @param X the design matrix
# #' @param y the trait vector
# #' @param tree the phylogenetic tree
# #'
# #' @details
# #' Initialize the stationary variance by taking the mad over the centered and
# #' normalized data. This assumes that the process as the tips is at the
# #' stationary state, so that all tips have variance gamma^2 = sigma^2 / (2 alpha).
# #' The data is centered using a (non phylogenetic) robust regression of y
# #' against X, and extracting the fitted values, that are taken to approximate the center.
# #'
# #' @return An initial gamma^2 parameter guess
# #'
# #' @keywords internal
# #'
# initAlphaParameter <- function(tree, trait, gamma_hat, tol = 0.1) {
#   ## Tree distance matrices
#   dist_phylo <- cophenetic.phylo(tree)
#   vcv_tree <- vcv(tree)
#   tree_heights <- diag(vcv_tree)
#   ## Keep only tol % closest tip paris
#   thershold <- quantile(dist_phylo[upper.tri(dist_phylo, diag = FALSE)], tol)
#   tip_pairs <- which(upper.tri(dist_phylo) & dist_phylo <= thershold, arr.ind = TRUE)
#   ## Normalize tip pairs differences
#   diff_tips <- abs(trait[tip_pairs[, 1]] - trait[tip_pairs[, 2]])
#   diff_tips <- diff_tips - disp_hat * (tree_heights[tip_pairs[, 1]] + tree_heights[tip_pairs[, 2]])
#   diff_tips <- diff_tips / (- 2 * disp_hat * apply(tip_pairs, 1, function(tt) vcv_tree[tt[1], tt[2]]))
#   ## Estimation
#   lambda_estim <- median(diff_tips)
#   lambda_estim <- abs(lambda_estim)
#   ## Bounds
#   max_lambda <- maxLambda(tree)
#   # if (lambda_estim < 0) return(0)
#   if (lambda_estim > max_lambda) return(max_lambda)
#   return(lambda_estim)
# }
#
# ##
# #' @title Initialization the selection strength alpha using robust estimation
# #'
# #' @description
# #' \code{init.alpha.estimation} fits (Y_i-Y_j)^2 ~ gamma^2(1-exp(-alpha*d_ij))
# #' for all couples of tips (i,j) that have the same mean, i.e than are not
# #' separated by a shift. Shifts are initialized thanks to a lasso
# #' (function \code{init.EM.lasso}).
# #'
# #' @details
# #' Function \code{robustbase::nlrob} is used for the robust fit.
# #'
# #' @param phylo a phylogenetic tree, class \code{\link[ape]{phylo}}.
# #' @param Y_data data at the tips.
# #' @param nbr_of_shifts : number of shifts wanted
# #' @param distances_phylo (matrix) : phylogenetic distance, result of function
# #' \code{compute_dist_phy}
# #' @param nbr_of_shifts number of shifts used in the EM algorithm
# #'
# #' @return params_init the list of initial parameters to be used, in the right
# #'  format.
# #'
# #' @keywords internal
# #'
# #10/07/14 - Initial release
# ##
# init.alpha.gamma.estimation <- function(phylo,
#                                         Y_data,
#                                         nbr_of_shifts,
#                                         times_shared,
#                                         distances_phylo,
#                                         T_tree,
#                                         subtree.list,
#                                         max_triplet_number,
#                                         alpha_known,
#                                         method.init.alpha.estimation,
#                                         tol_EM, h_tree,
#                                         miss,
#                                         masque_data,
#                                         independent, ...){
#   ## Initialize a vector with the group of each tip
#   tips_groups <- rep(0, length(phylo$tip.label))
#   names(tips_groups) <- phylo$tip.label
#   p <- nrow(Y_data)
#   ## Initialize shifts by a lasso without sigma
#   if (nbr_of_shifts > 0) {
#     lasso <- init.EM.lasso(phylo = phylo,
#                            Y_data = Y_data,
#                            process = "OU",
#                            nbr_of_shifts = nbr_of_shifts,
#                            use_sigma = FALSE,
#                            random.init = TRUE,
#                            stationary.root.init = TRUE,
#                            times_shared = times_shared,
#                            distances_phylo = distances_phylo,
#                            T_tree = T_tree,
#                            subtree.list = subtree.list,
#                            miss = miss,
#                            # impute_init_Rphylopars = FALSE,
#                            masque_data = masque_data,
#                            independent = independent,
#                            selection.strength.init = rep(1, p))
#     ## Roeorder phylo and trace edges
#     phy <- reorder(phylo, order = "cladewise")
#     edges_shifts <- correspondanceEdges(edges = lasso$shifts$edges,
#                                         from = phylo, to = phy)
#     ## Set groups of tips (one group = all the tips under a given shift)
#     Tr <- incidence.matrix(phy)
#     for (ed in order(edges_shifts)) { # Do it in order so that edges with higher numbers erase groups (edges closer from the tips)
#       ed_sh <- edges_shifts[ed]
#       tips_groups[phy$tip.label[Tr[, ed_sh]]] <- ed
#     }
#   } else {
#     edges_shifts <- NULL
#   }
#   ## For each group, take all the triplets of tips to estimate the covariance sigma_ij
#   cor_hat <- NULL # estimations from trilpets of pairs corelations
#   square_diff <- vector("list", p) # (Y_i-Y_j)^2
#   dists <- NULL # corresponding phylogenetic distances between pairs
#   hat_gam <- matrix(NA, nrow = length(edges_shifts)+1, ncol = p)
#   hat_gam_mad <- matrix(NA, nrow = length(edges_shifts)+1, ncol = p)
#   for (grp in 0:length(edges_shifts)) {
#     tips <- which(tips_groups==grp)
#     if (length(tips) > 1){
#       for (l in 1:p){
#         hat_gam[grp+1, l] <- var(na.omit(Y_data[l, tips]))
#         hat_gam_mad[grp+1, l] <- mad(Y_data[l, tips], na.rm = TRUE)^2
#         Z <- outer(Y_data[l, tips], Y_data[l, tips],
#                    function(x,y){x-y} )
#         square_diff[[l]] <- c(square_diff[[l]], (Z[upper.tri(Z)])^2)
#       }
#       Z <- distances_phylo[tips,tips]
#       dists <- c(dists, Z[upper.tri(Z)])
#     }
#   }
#   ## Estimation of gamma
#   gamma_0 <- matrix(NA, nrow = length(method.init.alpha.estimation) + 2, ncol = p)
#   rownames(gamma_0) <- c("var", "mad", method.init.alpha.estimation)
#   gamma_0["var", ] <- colMeans(hat_gam, na.rm = TRUE) # Simple variance
#   gamma_0["mad", ] <- robustbase::colMedians(hat_gam_mad, na.rm = TRUE) # MAD
#
#   ## Estimation of alpha
#   # Supress couple "too far away"
#   too_far <- (dists > h_tree)
#   dists <- dists[!too_far]
#   square_diff <- do.call(rbind, square_diff)
#   square_diff <- square_diff[, !too_far, drop = FALSE]
#   if (alpha_known) {
#     return(list(alpha_0 = init.alpha.gamma.default(alpha_known, ...)$alpha_0,
#                 gamma_0 = gamma_0[c("var", "mad")]))
#   } else {
#     alpha_0 <- matrix(NA, nrow = length(method.init.alpha.estimation), ncol = p)
#     rownames(alpha_0) <- method.init.alpha.estimation
#     for (method in method.init.alpha.estimation){
#       estimate.alpha  <- switch(method,
#                                 regression = estimate.alpha.regression,
#                                 regression.MM = estimate.alpha.regression.MM,
#                                 median = estimate.alpha.median)
#
#       for (l in 1:p){
#         mask <- !is.na(square_diff[l, ])
#         ag_0_try <- try(estimate.alpha(square_diff[l, mask],
#                                        dists[mask],
#                                        gamma_0["mad", l],
#                                        tol_EM, h_tree), silent = TRUE)
#
#         if (inherits(ag_0_try, "try-error")) {
#           message(paste0("Robust estimation of alpha by ", method, " failed."))
#           alpha_0[method, l] <- NA # init.alpha.gamma.default(alpha_known, ...)$alpha_0
#           gamma_0[method, l] <- NA
#         } else {
#           alpha_0[method, l] <- ag_0_try[["alpha_0"]]
#           gamma_0[method, l] <- ag_0_try[["gamma_0"]]
#         }
#       }
#     }
#     return(list(alpha_0 = alpha_0,
#                 gamma_0 = gamma_0))
#   }
# }
