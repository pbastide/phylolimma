#' @title Normalize RNASeq count data using gene lengths
#'
#' @description
#' This function normalizes a count matrix, using the matrix of length, and
#' an appropriate transformation.
#'
#' @param countMatrix The RNASeq count matrix. Rows and columns should be named.
#' @param lengthMatrix The associated length matrix. Should have the same dimensions as \code{countMatrix}, with the same names.
#' @param normalisationFactor Normalization factors to scale the raw library sizes, as computed e.g. by \code{\link[edgeR]{calcNormFactors}}.
#' @param lengthNormalization one of "none" (no correction), "TPM" (default) or "RPKM". See details.
#' @param dataTransformation one of "log2", "asin(sqrt)" or "sqrt." See details.
#'
#' @return A matrix of normalized and transformed counts, with the same dimensions as \code{countMatrix}.
#'
#' @details
#' The normalization procedures are:
#' \describe{
#' \item{\code{none}:}{No length normalization.}
#' \item{\code{TPM}:}{TODO}
#' \item{\code{RPKM}:}{TODO}
#' }
#'
#' @export
lengthNormalizeRNASeq <- function(countMatrix,
                                  lengthMatrix = NULL,
                                  normalisationFactor = NULL,
                                  lengthNormalization = c("TPM", "RPKM", "none"),
                                  dataTransformation = c("log2", "sqrt", "asin(sqrt)")) {

  ## Arguments
  lengthNormalization <- match.arg(lengthNormalization)
  dataTransformation <- match.arg(dataTransformation)

  ## Check matrices
  if (!is.matrix(countMatrix)) stop("'countMatrix' should be a matrix.")
  if (is.null(colnames(countMatrix))) stop("Column of count matrix should be named.")
  if (is.null(rownames(countMatrix))) stop("Rows of count matrix should be named.")

  if (lengthNormalization != "none") {
    if (!is.matrix(lengthMatrix)) stop("'lengthMatrix' should be a matrix.")
    if (any(dim(countMatrix) != dim(lengthMatrix))) stop("Count and length matrices should have the same dimension.")
    if (is.null(colnames(lengthMatrix))) stop("Column of length matrix should be named.")
    if (is.null(rownames(lengthMatrix))) stop("Rows of length matrix should be named.")
    if (any(rownames(countMatrix) != rownames(lengthMatrix))) stop("Count and length matrices should have the same row names.")
    if (any(colnames(countMatrix) != colnames(lengthMatrix))) stop("Count and length matrices should have the same column names.")
  }

  ## Check normalization factor
  if (is.null(normalisationFactor)) normalisationFactor <- 1.0
  if (!is.vector(normalisationFactor)) stop("'normalisationFactor' must be a vector.")
  if (length(normalisationFactor) > 1) {
    if (length(normalisationFactor) != ncol(countMatrix)) stop("'normalisationFactor' is a vector. Its length should be equal to the number of columns in 'countMatrix'.")
    if (is.null(names(normalisationFactor)) || any(names(normalisationFactor) != colnames(countMatrix))) stop("'normalisationFactor' is a vector. Its names should match the names of columns in 'countMatrix'.")
  }

  ## Normalization
  data.norm <- switch(lengthNormalization,
                      none = normalize_none(countMatrix, normalisationFactor, dataTransformation),
                      TPM = normalize_TPM(countMatrix, lengthMatrix, normalisationFactor, dataTransformation),
                      RPKM = normalize_RPKM(countMatrix, lengthMatrix, normalisationFactor, dataTransformation))

  ## Transformation
  data.trans <- switch(dataTransformation,
                       log2 = log2(data.norm),
                       "asin(sqrt)" = asin(sqrt(data.norm)),
                       sqrt = sqrt(data.norm))

  rownames(data.trans) <- rownames(countMatrix)
  colnames(data.trans) <- colnames(countMatrix)

  return(data.trans)
}

#' @title Normalize RNASeq count data
#'
#' @description
#' Apply standard CPM, with no length normalization.
#'
#' @inheritParams lengthNormalizeRNASeq
#'
#' @return A matrix of normalized count, with the same dimensions as \code{countMatrix}.
#'
#' @keywords internal
#'
normalize_none <- function(countMatrix, normalisationFactor, dataTransformation) {

  lib.size <- colSums(countMatrix) * normalisationFactor

  if (dataTransformation == "log2") {
    data.norm <- sweep(countMatrix + 0.5, 2, lib.size + 1, '/')
  } else {
    data.norm <- sweep(countMatrix, 2, lib.size, '/')
  }

  if (dataTransformation != "asin(sqrt)") data.norm <- data.norm * 1e6

  return(data.norm)
}

#' @title Normalize RNASeq count data
#'
#' @description
#' Apply TPM length normalization.
#'
#' @inheritParams lengthNormalizeRNASeq
#'
#' @return A matrix of normalized count, with the same dimensions as \code{countMatrix}.
#'
#' @keywords internal
#'
normalize_TPM <- function(countMatrix, lengthMatrix,
                          normalisationFactor, dataTransformation) {

  lib.size <- colSums(countMatrix / lengthMatrix) * normalisationFactor

  if (dataTransformation == "log2") {
    data.norm <- sweep((countMatrix + 0.5) / lengthMatrix, 2, lib.size + 1, '/')
  } else {
    data.norm <- sweep((countMatrix) / lengthMatrix, 2, lib.size, '/')
  }

  if (dataTransformation != "asin(sqrt)") data.norm <- data.norm * 1e6

  return(data.norm)
}

#' @title Normalize RNASeq count data
#'
#' @description
#' Apply RPKM length normalization.
#'
#' @inheritParams lengthNormalizeRNASeq
#'
#' @return A matrix of normalized count, with the same dimensions as \code{countMatrix}.
#'
#' @keywords internal
#'
normalize_RPKM <- function(countMatrix, lengthMatrix,
                           normalisationFactor, dataTransformation) {

  lib.size <- colSums(countMatrix) * normalisationFactor

  if (dataTransformation == "log2") {
    data.norm <- sweep((countMatrix + 0.5) / lengthMatrix, 2, lib.size + 1, '/')
  } else {
    data.norm <- sweep((countMatrix) / lengthMatrix, 2, lib.size, '/')
  }

  if (dataTransformation != "asin(sqrt)") data.norm <- data.norm * 1e9

  return(data.norm)
}
