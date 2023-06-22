test_that("phylolmFit - equivalencies", {
  set.seed(12891026)
  ## Tree
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)
  mat_tree <- ape::vcv(tree)
  ## data
  ngenes <- 20
  y_data <- t(phylolm::rTrait(ngenes, tree, model = "delta", parameters = list(delta = 0.1)))
  ## Design
  design <- matrix(1, nrow = ntips, ncol = 2)
  design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
  colnames(design) <- c("(Intercept)", "condition")
  rownames(design) <- tree$tip.label

  ## Test function
  test_lmFit_phylolm <- function(y_data, design, tree, model, measurement_error, REML) {
    ## Transform design and data and fit
    resLmFit <- phylolmFit(y_data, design = design, phy = tree,
                           model = model,
                           measurement_error = measurement_error,
                           use_consensus = FALSE,
                           REML = REML)

    ## phylolm fit
    bounds_alpha <- getBoundsSelectionStrength(tree)
    min_sig_err <- getMinError(tree)
    phylolm_fit <- function(y, design, phy, model, measurement_error, REML, ...) {
      data_phylolm <- as.data.frame(cbind(y, design))
      colnames(data_phylolm)[1] <- "expr"
      fplm <- suppressWarnings(phylolm::phylolm(expr ~ -1 + .,
                                                data = data_phylolm, phy = phy, model = model,
                                                measurement_error = measurement_error,
                                                upper.bound = list(alpha = bounds_alpha[2], lambda = getMaxLambda(min_sig_err)),
                                                lower.bound = list(alpha = bounds_alpha[1], sigma2_error = min_sig_err),
                                                REML = REML,
                                                ...))
      return(fplm)
    }
    fplm <- apply(y_data, 1, phylolm_fit, design = design, phy = tree, model = model, measurement_error = measurement_error, REML = REML)

    ## Comparison
    # coefficients
    phylolm_coef <- t(sapply(fplm, function(z) z$coefficients))
    expect_equivalent(phylolm_coef, resLmFit$coefficients, 1e-7)
    # df
    phylolm_df <- sapply(fplm, function(z) z$n - z$d)
    expect_equivalent(phylolm_df, resLmFit$df.residual)
    if (model == "BM" && !measurement_error) {
      # sigma
      phylolm_sigma <- sapply(fplm, function(z) sqrt(ntips/(ntips-2*(1-REML)) * z$sigma2))
      expect_equal(phylolm_sigma, resLmFit$sigma)
      # stdev
      phylolm_stdev <- t(sapply(fplm, function(z) sqrt(diag(z$vcov) / (ntips / (ntips-2*(1-REML)) * z$sigma2))))
      expect_equivalent(phylolm_stdev, resLmFit$stdev.unscaled)
    }
    # t.value
    phylolm_tvalue <- t(sapply(fplm, function(z) summary(z)$coefficients[, "t.value"]))
    phylolimma_tvalue <- resLmFit$coef / resLmFit$stdev.unscaled / resLmFit$sigma
    expect_equivalent(phylolm_tvalue, phylolimma_tvalue, 1e-7)
    # p.value
    phylolm_pvalue <- t(sapply(fplm, function(z) summary(z)$coefficients[, "p.value"]))
    phylolimma_pvalue <- 2 * pt(-abs(phylolimma_tvalue), df = resLmFit$df.residual)
    expect_equivalent(phylolm_pvalue, phylolimma_pvalue, 1e-7)
  }

  ## Tests
  model <- "BM"
  measurement_error <- FALSE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  measurement_error <- TRUE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)

  model <- "OUfixedRoot"
  measurement_error <- FALSE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  measurement_error <- TRUE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)

  model <- "lambda"
  measurement_error <- FALSE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)


  model <- "delta"
  measurement_error <- FALSE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  measurement_error <- TRUE
  REML <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)
  REML <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error, REML)

})
