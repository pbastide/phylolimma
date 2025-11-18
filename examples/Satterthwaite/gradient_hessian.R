#' @title Get Inverse Hessian
#'
#' @description
#' Compute the inverse Hessian
#'
#' @details
#' Code adapted from \code{lmerTest}, see
#' \url{https://github.com/runehaubo/lmerTestR/blob/35dc5885205d709cdc395b369b08ca2b7273cb78/R/lmer.R#L173}
#'
#' @param optpars parameter vector around which to compute the Hessian
#' @param fun function for which the Hessian neeeds to be computed
#' @param grad_trans gradient vector for transformed parameters
#' @param tol tolerence for calling zero eigenvalues
#'
#' @return The inverse hessian
#'
#' @keywords internal
#'
compute_inv_hessian <- function(optpars, fun, grad_trans, tol = 1e-8, ...) {
  # Compute Hessian:
  h <- numDeriv::hessian(func = fun, x = optpars, ...)
  # back transformation of parameters
  h <- t(grad_trans) %*% h %*% grad_trans
  # Eigen decompose the Hessian:
  eig_h <- eigen(h, symmetric=TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  zero <- evals > -tol & evals < tol
  if(sum(neg) > 0) { # negative eigenvalues
    eval_chr <- if(sum(neg) > 1) "eigenvalues" else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf("Model failed to converge with %d negative %s: %s",
                    sum(neg), eval_chr, evals_num), call.=FALSE)
  }
  # Note: we warn about negative AND zero eigenvalues:
  if(sum(zero) > 0) { # some eigenvalues are zero
    eval_chr <- if(sum(zero) > 1) "eigenvalues" else "eigenvalue"
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf("Model may not have converged with %d %s close to zero: %s",
                    sum(zero), eval_chr, evals_num))
  }
  # Compute vcov(varpar):
  pos <- eig_h$values > tol
  q <- sum(pos)
  # Using the Moore-Penrose generalized inverse for h:
  h_inv <- with(eig_h, {
    vectors[, pos, drop=FALSE] %*% diag(1/values[pos], nrow=q) %*%
      t(vectors[, pos, drop=FALSE]) })
  return(h_inv)
}
