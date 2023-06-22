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
                              use_consensus = FALSE,
                              ddf_method = "Samples")
  ## Fit
  resLmFit <- limma::lmFit(y_data, design = design)

  ## Test names and dimensions
  expect_equal(colnames(resPhyloLmFit$coefficients), colnames(resLmFit$coefficients))
  expect_equal(dim(resPhyloLmFit$coefficients), dim(resLmFit$coefficients))
  expect_equal(length(resPhyloLmFit$df.residual), length(resLmFit$df.residual))
  expect_equal(length(resPhyloLmFit$sigma), length(resLmFit$sigma))
  expect_equal(colnames(resPhyloLmFit$stdev.unscaled), colnames(resLmFit$stdev.unscaled))
  expect_equal(dim(resPhyloLmFit$stdev.unscaled), dim(resLmFit$stdev.unscaled))

  ## ddf
  expect_equal(resPhyloLmFit$df.residual, resLmFit$df.residual)
  expect_equal(getSpeciesNumber(tree), ntips)

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
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "BM",
                              measurement_error = TRUE,
                              lower.bound = list(sigma2_error = 0.1),
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.1 * resPhyloLmFit$sigma2_phy))

  ## Default lower bound
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "BM",
                              measurement_error = TRUE,
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= (.Machine$double.eps)^0.9 * resPhyloLmFit$sigma2_phy))

  ## Default lower bound
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "OUfixedRoot",
                              measurement_error = TRUE,
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= (.Machine$double.eps)^0.9 * resPhyloLmFit$sigma2_phy))

  ## Set upper bound and starting values
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "BM",
                              measurement_error = TRUE,
                              lower.bound = list(sigma2_error = 0.01),
                              upper.bound = list(sigma2_error = 0.5),
                              starting.value = list(sigma2_error = 0.05),
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.5 * resPhyloLmFit$sigma2_phy))

  ## Set upper bound and starting values
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "OUfixedRoot",
                              measurement_error = TRUE,
                              lower.bound = list(sigma2_error = 0.01),
                              upper.bound = list(sigma2_error = 0.5),
                              starting.value = list(sigma2_error = 0.05),
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.5 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))

  ## Set upper bound and starting values
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "OUfixedRoot",
                              measurement_error = TRUE,
                              lower.bound = list(alpha = 0.1, sigma2_error = 0.01),
                              upper.bound = list(alpha = 10, sigma2_error = 0.5),
                              starting.value = list(alpha = 5, sigma2_error = 0.05),
                              use_consensus = FALSE)

  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.5 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$optpar >= 0.1 - 1e-10))
  expect_true(all(resPhyloLmFit$optpar <= 10))

  ## Set upper bound and starting values
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = "OUfixedRoot",
                              measurement_error = TRUE,
                              lower.bound = list(sigma2_error = 0.01),
                              upper.bound = list(sigma2_error = 0.02),
                              starting.value = list(alpha = 1, sigma2_error = 0.015),
                              use_consensus = FALSE)

  tol <- 1e-16
  expect_true(all(resPhyloLmFit$sigma2_error >= 0.01 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar)))
  expect_true(all(resPhyloLmFit$sigma2_error <= 0.02 * resPhyloLmFit$sigma2_phy / (2 * resPhyloLmFit$optpar) + tol))

})

test_that("phylolmFit - mammals", {
  set.seed(12891026)
  tree <- ape::read.tree(text = "(((opossum.2:2e-12,opossum.1:0):2e-12,opossum.0:0):147.1,(armadillo.0:98.7,((((cow.2:2e-12,cow.1:0):2e-12,cow.0:0):85.4,(ferret.0:55.9,dog.0:55.9):29.5):10.8,((rabbit.0:89,(((rat.2:2e-12,rat.1:0):2e-12,rat.0:0):29.4,((((musCaroli.3:2e-12,musCaroli.2:0):2e-12,musCaroli.1:0):2e-12,musCaroli.0:0):4.1,(((((musMusculus.5:2e-12,musMusculus.4:0):2e-12,musMusculus.3:0):2e-12,musMusculus.2:0):2e-12,musMusculus.1:0):2e-12,musMusculus.0:0):4.1,((musSpretus.2:2e-12,musSpretus.1:0):2e-12,musSpretus.0:0):4.1):25.3):59.6):2.9,((marmoset.1:2e-12,marmoset.0:0):52.4,(((((rhesus.4:2e-12,rhesus.3:0):2e-12,rhesus.2:0):2e-12,rhesus.1:0):2e-12,rhesus.0:0):35.3,((orangutan.1:2e-12,orangutan.0:0):18.4,((gorilla.1:2e-12,gorilla.0:0):11.5,((human.1:2e-12,human.0:0):8.8,((bonobo.1:2e-12,bonobo.0:0):4.5,(chimp.1:2e-12,chimp.0:0):4.5):4.3):2.7):6.9):16.9):17.1):39.5):4.3):2.5):48.4);")
  ntips <- length(tree$tip.label)
  h_tree <- max(ape::vcv(tree))
  ## data
  ngenes <- 20
  y_data <- t(phylolm::rTrait(ngenes, tree, model = "lambda", parameters = list(lambda = 1)))
  ## Design
  design <- matrix(1, nrow = ntips, ncol = 2)
  design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
  colnames(design) <- c("(Intercept)", "condition")
  rownames(design) <- tree$tip.label

  min_sig_err <- getMinError(tree)

  pp <- phylolmFit(y_data, design = design, phy = tree,
                   model = "BM",
                   lower.bound = list(sigma2_error = 0.001),
                   measurement_error = TRUE,
                   use_consensus = FALSE,
                   ddf_method = "Samples")
  expect_true(all(pp$sigma2_error >= 0.001 * pp$sigma2_phy))

  expect_equal(getSpeciesNumber(tree), length(unique(sub("\\.[0-9]", "", tree$tip.label))))

  pp <- phylolmFit(y_data, design = design, phy = tree,
                   model = "BM",
                   measurement_error = TRUE,
                   use_consensus = FALSE,
                   ddf_method = "Samples")
  expect_true(all(pp$sigma2_error >= min_sig_err * pp$sigma2_phy))

  pp <- phylolmFit(y_data, design = design, phy = tree,
                   model = "OUfixedRoot",
                   lower.bound = list(alpha = 100),
                   upper.bound = list(alpha = 101),
                   starting.value = list(alpha = 100.5),
                   measurement_error = TRUE,
                   use_consensus = FALSE)

  expect_true(all(pp$sigma2_error >= min_sig_err * pp$sigma2_phy / (2 * pp$optpar)))

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
