test_that("transformation is correct", {
  set.seed(12891026)
  ## Tree
  ntips <- 50
  tree <- ape::rphylo(ntips, 0.1, 0)
  mat_tree <- ape::vcv(tree)
  ## Chol
  C_tree <- get_chol_tree(NULL, NULL, tree, "BM", FALSE)$C_tree
  expect_equal(t(C_tree) %*% C_tree, mat_tree)
  ## data
  ngenes <- 100
  y_data <- matrix(rnorm(ngenes * ntips, 0, 1), ncol = ntips)
  design <- matrix(1, nrow = ntips, ncol = 2)
  design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
  ## Trans
  design_trans <- transform_design_tree(C_tree, design)
  expect_equal(dim(design_trans), dim(design))
  y_trans <- transform_data_tree(C_tree, y_data)
  expect_equal(dim(t(y_trans)), dim(y_data))
})

test_that("phylolmFit - BM", {
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
  ## Model
  model <- "BM"
  measurement_error <- TRUE
  ## Fit Phylo
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = model,
                              measurement_error = measurement_error,
                              use_consensus = FALSE)
  ## Fit
  resLmFit <- limma::lmFit(y_data, design = design)

  ## Test names and dimensions
  expect_equal(colnames(resPhyloLmFit$coefficients), colnames(resLmFit$coefficients))
  expect_equal(dim(resPhyloLmFit$coefficients), dim(resLmFit$coefficients))
  expect_equal(length(resPhyloLmFit$df.residual), length(resLmFit$df.residual))
  expect_equal(length(resPhyloLmFit$sigma), length(resLmFit$sigma))
  expect_equal(colnames(resPhyloLmFit$stdev.unscaled), colnames(resLmFit$stdev.unscaled))
  expect_equal(dim(resPhyloLmFit$stdev.unscaled), dim(resLmFit$stdev.unscaled))

  ## ebayes
  fitphy <- eBayes(resPhyloLmFit)
  fitlimma <- limma::eBayes(resLmFit)

  ## Test names and dimensions
  expect_equal(colnames(fitphy$p.value), colnames(fitlimma$p.value))
  expect_equal(dim(fitphy$p.value), dim(fitlimma$p.value))

  ## Other functions
  expect_error(treat(resPhyloLmFit), "is not supported for an object of class `PhyloMArrayLM`.")
  expect_error(decideTests(resPhyloLmFit), "is not supported for an object of class `PhyloMArrayLM`.")
  expect_error(classifyTestsF(resPhyloLmFit), "is not supported for an object of class `PhyloMArrayLM`.")

})

test_that("phylolmFit - bounds", {
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

  ## Set lower bound
  resPhyloLmFit <- expect_warning(
    phylolmFit(y_data, design = design, phy = tree,
               model = "BM",
               measurement_error = TRUE,
               lower.bound = list(sigma2_error = 0.1),
               use_consensus = FALSE),
    "the estimation of sigma2_error matches the upper/lower bound for this parameter.")

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.1 * resPhyloLmFit$sigma2_phy))

  ## Default lower bound
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "BM",
                              measurement_error = TRUE,
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= (.Machine$double.eps)^0.9 * resPhyloLmFit$sigma2_phy))

  ## Default lower bound
  resPhyloLmFit <- expect_warning(
    phylolmFit(y_data, design = design, phy = tree,
               model = "OUfixedRoot",
               measurement_error = TRUE,
               use_consensus = FALSE),
    "the estimation of sigma2_error matches the upper/lower bound for this parameter.")

  expect_true(all(resPhyloLmFit$sigma2_error >= (.Machine$double.eps)^0.9 * resPhyloLmFit$sigma2_phy))

  ## Set upper bound and starting values
  resPhyloLmFit <- expect_warning(
    phylolmFit(y_data, design = design, phy = tree,
               model = "BM",
               measurement_error = TRUE,
               lower.bound = list(sigma2_error = 0.01),
               upper.bound = list(sigma2_error = 0.5),
               starting.value = list(sigma2_error = 0.05),
               use_consensus = FALSE),
    "the estimation of sigma2_error matches the upper/lower bound for this parameter.")

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.5 * resPhyloLmFit$sigma2_phy))

  ## Set upper bound and starting values
  resPhyloLmFit <- expect_warning(
    phylolmFit(y_data, design = design, phy = tree,
               model = "OUfixedRoot",
               measurement_error = TRUE,
               lower.bound = list(sigma2_error = 0.01),
               upper.bound = list(sigma2_error = 0.5),
               starting.value = list(sigma2_error = 0.05),
               use_consensus = FALSE),
    "the estimation of sigma2_error matches the upper/lower bound for this parameter.")

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.5 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))

  ## Set upper bound and starting values
  expect_warning(resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                                             model = "OUfixedRoot",
                                             measurement_error = TRUE,
                                             lower.bound = list(alpha = 0.1, sigma2_error = 0.01),
                                             upper.bound = list(alpha = 10, sigma2_error = 0.5),
                                             starting.value = list(alpha = 5, sigma2_error = 0.05),
                                             use_consensus = FALSE))

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.5 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$optpar >= 0.1))
  expect_true(all(resPhyloLmFit$optpar <= 10))

  ## Set upper bound and starting values
  resPhyloLmFit <- expect_warning(
    phylolmFit(y_data, design = design, phy = tree,
               model = "OUfixedRoot",
               measurement_error = TRUE,
               lower.bound = list(sigma2_error = 0.01),
               upper.bound = list(sigma2_error = 0.02),
               starting.value = list(alpha = 1, sigma2_error = 0.015),
               use_consensus = FALSE),
    "the estimation of sigma2_error matches the upper/lower bound for this parameter.")

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.02 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))

})

# test_that("phylolmFit - BM - error", {
#   set.seed(12891026)
#   ## Tree
#   ntips <- 50
#   tree <- ape::rphylo(ntips, 0.1, 0)
#   ## data
#   ngenes <- 100
#   y_data <- matrix(rnorm(ngenes * ntips, 0, 1), ncol = ntips)
#   design <- matrix(1, nrow = ntips, ncol = 2)
#   design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
#   colnames(y_data) <- rownames(design) <- tree$tip.label
#   ## Fit
#   fit <- phylolmFit(y_data, design = design, phy = tree,
#                     model = "BM", measurement_error = TRUE,
#                     ndups = 1, spacing = 1, block = NULL, weights = NULL, method = "ls")
#   ## ebayes
#   fitb <- limma::eBayes(fit)
# })
