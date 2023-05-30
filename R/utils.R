#' @title Get Bounds on alpha
#'
#' @description
#' Find reasonable bounds on the \code{alpha} parameter.
#'
#' @param phy a phylogenetic tree.
#' @param relative_half_life_min optional minimal half life relative to tree height
#' @param relative_half_life_max optional maximal half life relative to tree height
#'
#' @details
#' This functions tries to find reasonable bounds on the \eqn{\alpha} parameter
#' of an OU process by using the scaled phylogenetic half-life \eqn{t_{1/2} = \log(2) / \alpha / h},
#' where \eqn{h} is the total height of the tree.
#' If \eqn{t_{1/2} = D)}, it means that the trait will need a time \eqn{D \times h}
#' to cover half the distance to the optimal value (see Hansen, 1997).
#' Small values of \eqn{D} means high selection pressure (large \eqn{\alpha}),
#' while large values of \eqn{D} means low selection pressure (small \eqn{\alpha}).
#'
#' The default maximum value for \eqn{D} is \code{relative_half_life_max = 10000}
#' (selection is week and the process looks like a BM).
#'
#' The default minimum value for \eqn{D} is \code{relative_half_life_min = 0.0001}
#' (selection is strong and tips are only weakly correlated).
#'
#' The function makes sure that the maximum \eqn{\alpha} value associated with
#' \code{relative_half_life_min} does not lead to underflow errors when rescaling
#' the tree.
#'
#' @return A vector with lower and upper bounds for alpha
#'
#' @references
#' Hansen, T. F. (1997). Stabilizing Selection and the Comparative Analysis of Adaptation. Evolution, 51(5) :1341.
#'
#' @keywords internal
#'
getBoundsSelectionStrength <- function(phy,
                                       relative_half_life_min = 0.0001,
                                       relative_half_life_max = 10000) {
  ## alpha min
  h_tree <- tree_height(phy)
  alpha_min <- log(2) / (relative_half_life_max * h_tree)
  ## alpha max
  alpha_max <- log(2) / (relative_half_life_min * h_tree)
  alpha_max_machine <- log(.Machine$double.xmax^0.975) / (2 * h_tree)
  if (alpha_max > alpha_max_machine) alpha_max <- alpha_max_machine
  return(c(alpha_min, alpha_max))
}

#' @title Get Min on sigma2_error
#'
#' @description
#' Find reasonable minimum on the \code{sigma2_error} parameter.
#'
#' @param phy a phylogenetic tree.
#' @param tol the numerical tolerance
#'
#' @details
#' The minimum value must be high enough so that it can be numerically
#' distinguished from zero.
#' Default to \eqn{tol * h}, where \eqn{h} is the total height of the tree.
#' If an OU is used, then this value is updated to match the transformed tree height.
#'
#' @return The minimum value for sigma2_error
#'
#' @keywords internal
#'
getMinError <- function(phy,
                        tol = (.Machine$double.eps)^0.5) {
  h_tree <- tree_height(phy)
  return(tol * h_tree)
}

#' @title Get Max on lambda
#'
#' @description
#' Find reasonable maximum on the \code{lambda} parameter.
#'
#' @param min_error inimum error value
#'
#' @return The maximum value for lambda
#'
#' @keywords internal
#'
getMaxLambda <- function(min_error) {
  return(1 / (1 + min_error))
}

#' @title Get lower bounds on parameters
#'
#' @description
#' Get the lower bounds on the parameters
#'
#' @param ... user specified parameters
#'
#' @return list of lower bounds on parameters
#'
#' @keywords internal
#'
get_lower_bounds <- function(bounds_alpha, min_sigma2_error, ...) {
  dot_args <- dots(...)
  lower_bounds <- list()
  if ("lower.bound" %in% names(dot_args)) lower_bounds <- eval(dot_args$lower.bound)
  # If no bounds specified on sigma2_error, set it.
  if (!("sigma2_error" %in% names(lower_bounds))) lower_bounds$sigma2_error = min_sigma2_error
  # If no bounds specified on alpha, set it.
  if (!("alpha" %in% names(lower_bounds))) lower_bounds$alpha = bounds_alpha[1]
  return(as.list(lower_bounds))
}

#' @title Get upper bounds on parameters
#'
#' @description
#' Get the upper bounds on the parameters
#'
#' @param ... user specified parameters
#'
#' @return list of upper bounds on parameters
#'
#' @keywords internal
#'
get_upper_bounds <- function(bounds_alpha, min_sigma2_error, ...) {
  dot_args <- dots(...)
  upper_bounds <- list()
  if ("upper.bound" %in% names(dot_args)) upper_bounds <- eval(dot_args$upper.bound)
  # If no bounds specified on alpha, set it.
  if (!("alpha" %in% names(upper_bounds))) upper_bounds$alpha = bounds_alpha[2]
  # If no bounds specified on lambda, set it.
  if (!("lambda" %in% names(upper_bounds))) upper_bounds$lambda = getMaxLambda(min_sigma2_error)
  return(as.list(upper_bounds))
}

#' @title Get starting values on parameters
#'
#' @description
#' Get the starting values of the parameters
#'
#' @param ... user specified parameters
#'
#' @return list of starting values on parameters
#'
#' @keywords internal
#'
get_starting_values <- function(bounds_alpha, ...) {
  dot_args <- dots(...)
  starting_values <- list()
  if ("starting.value" %in% names(dot_args)) starting_values <- eval(dot_args$starting.value)
  # If no bounds specified on alpha, set it.
  # if (!("alpha" %in% names(starting_values))) starting_values$alpha = mean(bounds_alpha)
  ## TODO init alpha from cheries ?
  return(as.list(starting_values))
}

#' @title Get all paramters but lower and upper bounds and starting value
#'
#' @description
#' Get all parameters but \code{lower.bounds}, \code{upper.bound} and \code{starting.value}
#'
#' @param ... user specified parameters
#'
#' @return list of parameters
#'
#' @keywords internal
#'
get_dots_args <- function(...) {
  dot_args <- dots(...)
  if ("lower.bound" %in% names(dot_args)) {
    dot_args <- dot_args[names(dot_args) != "lower.bound"]
  }
  if ("upper.bound" %in% names(dot_args)) {
    dot_args <- dot_args[names(dot_args) != "upper.bound"]
  }
  if ("starting.value" %in% names(dot_args)) {
    dot_args <- dot_args[names(dot_args) != "starting.value"]
  }
  return(dot_args)
}

#' @title Get maximum tree height
#'
#' @param tree a phylogenetic tree
#'
#' @return maximum tree height
#'
#' @keywords internal
#'
tree_height <- function(tree) {
  return(max(ape::node.depth.edgelength(tree)))
}

#' @title Get lambda value
#'
#' @description
#' Compute the lambda value equivalent to the sigma2 plus sigma2_error model.
#'
#'
#' @param sigma2 phylogenetic variance
#' @param sigma2_error measurement error variance
#' @param h_tree height of the tree
#'
#' @return lambda_error
#'
#' @keywords internal
#'
get_lambda_error <- function(sigma2, sigma2_error, h_tree) {
  return(sigma2 * h_tree / (sigma2_error + sigma2 * h_tree))
}
