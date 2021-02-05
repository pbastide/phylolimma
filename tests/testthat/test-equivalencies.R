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
  test_lmFit_phylolm <- function(y_data, design, tree, model, measurement_error) {
    ## Transform design and data and fit
    if (model != "delta") {
      resLmFit <- phylolmFit(y_data, design = design, phy = tree,
                             model = model,
                             measurement_error = measurement_error)
    } else {
      expect_warning(resLmFit <- phylolmFit(y_data, design = design, phy = tree,
                                            model = model,
                                            measurement_error = measurement_error),
                     "the estimation of delta matches the upper/lower bound for this parameter.")
    }

    ## phylolm fit
    phylolm_fit <- function(y, design, phy, model, measurement_error, ...) {
      data_phylolm <- as.data.frame(cbind(y, design))
      colnames(data_phylolm)[1] <- "expr"
      fplm <- phylolm::phylolm(expr ~ -1 + ., data = data_phylolm, phy = phy, model = model, measurement_error = measurement_error, ...)
      return(fplm)
    }
    if (model != "delta") {
      fplm <- apply(y_data, 1, phylolm_fit, design = design, phy = tree, model = model, measurement_error = measurement_error)
    } else {
      expect_warning(fplm <- apply(y_data, 1, phylolm_fit, design = design, phy = tree, model = model, measurement_error = measurement_error),
                     "the estimation of delta matches the upper/lower bound for this parameter.")
    }

    ## Comparison
    # coefficients
    phylolm_coef <- t(sapply(fplm, function(z) z$coefficients))
    expect_equivalent(phylolm_coef, resLmFit$coefficients)
    # df
    phylolm_df <- sapply(fplm, function(z) z$n - z$d)
    expect_equivalent(phylolm_df, resLmFit$df.residual)
    if (model == "BM" && !measurement_error) {
      # sigma
      phylolm_sigma <- sapply(fplm, function(z) sqrt(ntips/(ntips-2) * z$sigma2))
      expect_equal(phylolm_sigma, resLmFit$sigma)
      # stdev
      phylolm_stdev <- t(sapply(fplm, function(z) sqrt(diag(z$vcov) / (ntips / (ntips-2) * z$sigma2))))
      expect_equivalent(phylolm_stdev, resLmFit$stdev.unscaled)
    }
    # t.value
    phylolm_tvalue <- t(sapply(fplm, function(z) summary(z)$coefficients[, "t.value"]))
    phylolimma_tvalue <- resLmFit$coef / resLmFit$stdev.unscaled / resLmFit$sigma
    expect_equivalent(phylolm_tvalue, phylolimma_tvalue)
    # p.value
    phylolm_pvalue <- t(sapply(fplm, function(z) summary(z)$coefficients[, "p.value"]))
    phylolimma_pvalue <- 2 * pt(-abs(phylolimma_tvalue), df = resLmFit$df.residual)
    expect_equivalent(phylolm_pvalue, phylolimma_pvalue)
  }

  ## Tests
  model <- "BM"
  measurement_error <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)
  measurement_error <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)

  model <- "OUfixedRoot"
  measurement_error <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)
  measurement_error <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)

  model <- "lambda"
  measurement_error <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)


  model <- "delta"
  measurement_error <- FALSE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)
  measurement_error <- TRUE
  test_lmFit_phylolm(y_data, design, tree, model, measurement_error)

})
