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


#' @title Normalize RNASeq count data from tximport
#'
#' @description
#' Normalize RNASeq data from a tximport object using DESeq2.
#' This function handles length-scaled TPM counts from Salmon/kallisto via tximport,
#' applying DESeq2's size factor estimation and variance stabilizing transformations.
#'
#' @param txi A tximport object containing counts, abundance, and length matrices.
#' @param colData A data.frame with sample information. Row names must match column names of \code{txi$counts}.
#' @param design A formula specifying the design for DESeq2 (e.g., \code{~ condition}).
#' @param removeBatch Optional column name in \code{colData} specifying a batch variable to remove using \code{\link[limma]{removeBatchEffect}}.
#' @param dataTransformation One of "log2", "vst", "sqrt", or "asin(sqrt)". See details.
#'
#' @return A matrix of normalized and transformed counts.
#'
#' @details
#' This function uses DESeq2's \code{DESeqDataSetFromTximport} which properly handles
#' the length-scaled counts from tximport. Size factors are estimated using DESeq2's
#' median-of-ratios method.
#'
#' If \code{removeBatch} is specified, batch effects are removed using limma's
#' \code{removeBatchEffect} while preserving the design effects.
#'
#' The data transformations are:
#' \describe{
#' \item{\code{log2}:}{Log2 transformation with pseudo-count of 0.5.}
#' \item{\code{vst}:}{Variance stabilizing transformation from DESeq2.}
#' \item{\code{sqrt}:}{Square root transformation with pseudo-count of 0.5.}
#' \item{\code{asin(sqrt)}:}{Arcsine square root transformation with pseudo-count of 0.5.}
#' }
#'
#' @references
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data}
#'
#' @export
normalizeFromTximport <- function(txi,
                                  colData,
                                  design,
                                  removeBatch        = NULL,
                                  dataTransformation = c("log2", "vst", "sqrt", "asin(sqrt)")) {

  ## Arguments
  dataTransformation <- match.arg(dataTransformation)

  ## Check txi
 if (!is.list(txi) || !all(c("counts", "abundance", "length") %in% names(txi))) {
    stop("'txi' must be a tximport object with 'counts', 'abundance', and 'length' elements.")
  }

  ## Check colData
  if (!is.data.frame(colData)) stop("'colData' must be a data.frame.")
  if (!all(colnames(txi$counts) %in% rownames(colData))) {
    stop("Row names of 'colData' must match column names of 'txi$counts'.")
  }

  ## Check removeBatch
  if (!is.null(removeBatch) && !(removeBatch %in% colnames(colData))) {
    stop("'removeBatch' column '", removeBatch, "' not found in 'colData'.")
  }

  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = colData, design = design)
  dds <- DESeq2::estimateSizeFactors(dds)

  data.trans <- applyTransformation(dds, dataTransformation)

  if (!is.null(removeBatch)) {
    mm        <- model.matrix(design, colData)
    batch_col <- which(colnames(mm) == removeBatch)
    ## Keep the design matrix with all experimental factors other than the batch effects
    if (length(batch_col) > 0) mm <- mm[, -batch_col, drop = FALSE]
    data.trans <- limma::removeBatchEffect(data.trans, batch = colData[[removeBatch]], design = mm)
  }

  return(data.trans)
}


#' @title Normalize RNASeq count data from a count matrix
#'
#' @description
#' Normalize RNASeq data from a count matrix using DESeq2.
#' Optionally accounts for gene length using a length matrix.
#'
#' @param countMatrix The RNASeq count matrix. Rows are genes, columns are samples. Must be named.
#' @param colData A data.frame with sample information. Row names must match column names of \code{countMatrix}.
#' @param design A formula specifying the design for DESeq2 (e.g., \code{~ condition}).
#' @param removeBatch Optional column name in \code{colData} specifying a batch variable to remove using \code{\link[limma]{removeBatchEffect}}.
#' @param lengthMatrix Optional length matrix for length normalization. Should have the same dimensions and names as \code{countMatrix}.
#' @param dataTransformation One of "log2", "vst", "sqrt", or "asin(sqrt)". See details.
#'
#' @return A matrix of normalized and transformed counts.
#'
#' @details
#' When \code{lengthMatrix} is provided, it is normalized by dividing each row by its
#' geometric mean, then passed to DESeq2's \code{estimateSizeFactors} as \code{normMatrix}.
#' This accounts for gene length bias while using DESeq2's median-of-ratios normalization.
#'
#' If \code{removeBatch} is specified, batch effects are removed using limma's
#' \code{removeBatchEffect} while preserving the design effects.
#'
#' The data transformations are:
#' \describe{
#' \item{\code{log2}:}{Log2 transformation with pseudo-count of 0.5.}
#' \item{\code{vst}:}{Variance stabilizing transformation from DESeq2.}
#' \item{\code{sqrt}:}{Square root transformation with pseudo-count of 0.5.}
#' \item{\code{asin(sqrt)}:}{Arcsine square root transformation with pseudo-count of 0.5.}
#' }
#'
#' @export
normalizeFromMatrix <- function(countMatrix,
                                colData,
                                design,
                                removeBatch        = NULL,
                                lengthMatrix       = NULL,
                                dataTransformation = c("log2", "vst", "sqrt", "asin(sqrt)")) {

  ## Arguments
  dataTransformation <- match.arg(dataTransformation)

  ## Check matrices
  if (!is.matrix(countMatrix))        stop("'countMatrix' should be a matrix.")
  if (is.null(colnames(countMatrix))) stop("Column of count matrix should be named.")
  if (is.null(rownames(countMatrix))) stop("Rows of count matrix should be named.")

  ## Check colData
  if (!is.data.frame(colData)) stop("'colData' must be a data.frame.")
  if (!all(colnames(countMatrix) %in% rownames(colData))) {
    stop("Row names of 'colData' must match column names of 'countMatrix'.")
  }

  ## Check removeBatch
  if (!is.null(removeBatch) && !(removeBatch %in% colnames(colData))) {
    stop("'removeBatch' column '", removeBatch, "' not found in 'colData'.")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countMatrix, colData = colData, design = design)

  if (!is.null(lengthMatrix)) {
    if (!is.matrix(lengthMatrix))                             stop("'lengthMatrix' should be a matrix.")
    if (any(dim(countMatrix) != dim(lengthMatrix)))           stop("Count and length matrices should have the same dimension.")
    if (is.null(colnames(lengthMatrix)))                      stop("Column of length matrix should be named.")
    if (is.null(rownames(lengthMatrix)))                      stop("Rows of length matrix should be named.")
    if (any(rownames(countMatrix) != rownames(lengthMatrix))) stop("Count and length matrices should have the same row names.")
    if (any(colnames(countMatrix) != colnames(lengthMatrix))) stop("Count and length matrices should have the same column names.")

    # DESeq2 models counts on the log scale with an intercept. Since gene length
    # affects counts multiplicatively, we center the length matrix by dividing
    # each row by its geometric mean (centering at zero on the log scale).
    lengthMatrix <- lengthMatrix / exp(rowMeans(log(lengthMatrix)))

    dds <- DESeq2::estimateSizeFactors(dds, normMatrix = lengthMatrix)
  } else {
    dds <- DESeq2::estimateSizeFactors(dds)
  }

  data.trans <- applyTransformation(dds, dataTransformation)

  if (!is.null(removeBatch)) {
    mm        <- model.matrix(design, colData)
    batch_col <- which(colnames(mm) == removeBatch)
    ## Keep the design matrix with all experimental factors other than the batch effects
    if (length(batch_col) > 0) mm <- mm[, -batch_col, drop = FALSE]
    data.trans <- limma::removeBatchEffect(data.trans, batch = colData[[removeBatch]], design = mm)
  }

  return(data.trans)
}

#' @title Apply data transformation to DESeq2 object

#' @description
#' Apply a variance-stabilizing or other transformation to normalized counts.
#'
#' @param dds A DESeqDataSet object with estimated size factors.
#' @param dataTransformation One of "log2", "vst", "sqrt", or "asin(sqrt)".
#'
#' @return A matrix of transformed counts.
#'
#' @keywords internal
applyTransformation <- function(dds, dataTransformation) {
  if (dataTransformation == "vst") {
    data.trans <- SummarizedExperiment::assay(DESeq2::vst(dds))
  } else {
    data.norm  <- DESeq2::counts(dds, normalized = TRUE)
    data.trans <- switch(dataTransformation,
         log2         = log2(data.norm      + 0.5),
         "asin(sqrt)" = asin(sqrt(data.norm + 0.5)),
         sqrt         = sqrt(data.norm      + 0.5))
  }

  return(data.trans)
}
