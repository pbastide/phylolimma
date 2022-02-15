#' @title Transform RNA-Seq Data Ready for Linear Modelling
#'
#' @description
#' Normalize and transform count data,
#' estimate the mean-variance relationship
#' and use this to compute appropriate observation-level weights.
#' The data are then ready for phylogenetic linear modeling.
#'
#' @inheritParams lengthNormalizeRNASeq
#' @inheritParams limma::voom
#' @inheritParams phylolmFit
#' @param ... further parameters to be passed to \code{\link{phylolmFit}}.
#'
#' @return A matrix of normalized and transformed counts, with the same dimensions as \code{counts}.
#'
#' @details
#' See details of the data normalization and transformation in function
#' \code{\link{lengthNormalizeRNASeq}}.
#'
#' @export
phyloVoom <- function(counts,
                      # >>>>>>>>>
                      lengthMatrix = NULL,
                      normalisationFactor = NULL,
                      lengthNormalization = c("TPM", "RPKM", "none"),
                      dataTransformation = "log2",
                      # <<<<<<<<<
                      design = NULL,
                      # >>>>>>>>>
                      phy,
                      model = c("BM", "lambda", "OUfixedRoot", "delta"),
                      measurement_error = TRUE,
                      # <<<<<<<<<
                      weights = NULL, span = 0.5, plot = FALSE, save.plot = FALSE,
                      ...)
  #	Linear modelling of count data with mean-variance modelling at the observation level.
  #	Creates an EList object for entry to lmFit() etc in the limma pipeline.
  #	Gordon Smyth and Charity Law
  #	Created 22 June 2011.  Last modified 23 January 2020.
  # >>>>>>>>>
  # Modified to include sequence length normalization and phylolm fit
  # Changes taged by ">>>>>>>>> ...  <<<<<<<<<" flags
  # Paul Bastide
  # <<<<<<<<<
{
  out <- list()

  #	Extract counts from known data objects
  if(is(counts,"DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
    if(is.null(normalisationFactor)) normalisationFactor <- counts$samples$norm.factors
    counts <- counts$counts
  } else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
    if(isExpressionSet) {
      if(length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
      if(length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    } else {
      counts <- as.matrix(counts)
    }
  }

  #	Check counts
  n <- nrow(counts)
  if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")
  m <- min(counts)
  if(is.na(m)) stop("NA counts not allowed")
  if(m < 0) stop("Negative counts now allowed")

  #	Check design
  if(is.null(design)) {
    design <- matrix(1,ncol(counts),1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }

  #>>>>>>>>> Begin Change
  #	Check lib.size
  if (is.null(normalisationFactor)) normalisationFactor <- 1.0
  lib.size <- colSums(counts) * normalisationFactor

  # Check model
  if (!measurement_error) stop("'phyloVoom' can only be applied on models with 'measurement_error=TRUE'.")
  if (dataTransformation != "log2") stop("'phyloVoom' can only be applied for log2 transformed models.")

  #	Fit linear model to normalized data
  y <- lengthNormalizeRNASeq(counts,
                             lengthMatrix = lengthMatrix,
                             normalisationFactor = normalisationFactor,
                             lengthNormalization = lengthNormalization,
                             dataTransformation = dataTransformation)
  # browser()
  fit <- phylolmFit(y, design, phy = phy, model = model, measurement_error = measurement_error, weights = weights, ...)
  if(is.null(fit$AmeanRaw)) fit$AmeanRaw <- rowMeans(y,na.rm=TRUE)
  #<<<<<<<<< End Change

  #	If no replication found, set all weight to 1
  NWithReps <- sum(fit$df.residual > 0L)
  if(NWithReps < 2L) {
    if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
    if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    out$design <- design
    if(is.null(out$targets))
      out$targets <- data.frame(lib.size=lib.size)
    else
      out$targets$lib.size <- lib.size
    return(new("EList",out))
  }
  #	Fit lowess trend to sqrt-standard-deviations by log-count-size
  #>>>>>>>>> Begin Change
  # browser()
  resx <- invert_log2_normalisation(fit$AmeanRaw, counts, lengthMatrix, normalisationFactor, lengthNormalization)
  sx <- resx$sx
  sy <- (fit$sigma2_error)^(1/4)
  #<<<<<<<<< End Change
  allzero <- rowSums(counts)==0
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx,sy,f=span)
  if(plot) {
    #>>>>>>>>> Begin Change
    plot(sx,sy,xlab=resx$name_sx,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
    #<<<<<<<<< End Change
    title("phyloVoom: Mean-variance trend")
    lines(l,col="red")
  }

  #	Make interpolating rule
  #	Special treatment of zero counts is now removed;
  #	instead zero counts get same variance as smallest gene average.
  #	l$x <- c(0.5^0.25, l$x)
  #	l$x <- c(log2(0.5), l$x)
  #	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
  #	var0 <- max(var0,1e-6)
  #	l$y <- c(var0, l$y)
  f <- approxfun(l, rule=2, ties=list("ordered",mean))

  #	Find individual quarter-root fitted counts
  #>>>>>>>>> Begin Change
  fitted.values <- predict_phylolmFit(fit)
  fitted.logcount <- invert_log2_fitted(fitted.values, counts, lengthMatrix, normalisationFactor, lengthNormalization)
  #<<<<<<<<< End Change

  #	Apply trend to individual observations
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  colnames(w) <- colnames(counts)

  #	Output
  out$E <- y
  out$weights <- w
  out$design <- design
  if(is.null(out$targets))
    out$targets <- data.frame(lib.size=lib.size)
  else
    out$targets$lib.size <- lib.size
  if(save.plot) {
    #>>>>>>>>> Begin Change
    out$voom.xy <- list(x=sx,y=sy,xlab=resx$name_sx,ylab="Sqrt( standard deviation of errors)")
    #<<<<<<<<< End Change
    out$voom.line <- l
  }

  new("EList",out)
}
