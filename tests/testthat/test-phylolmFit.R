test_that("transformation is correct", {
  set.seed(12891026)
  ## Tree
  ntips <- 50
  tree <- ape::rphylo(ntips, 0.1, 0)
  mat_tree <- ape::vcv(tree)
  ## Chol
  C_tree <- get_chol_tree(NULL, NULL, tree, "BM", FALSE)
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

# test_that("phylolmFit - BM", {
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
#                     model = "BM",
#                     ndups = 1, spacing = 1, block = NULL, weights = NULL, method = "ls")
#   ## ebayes
#   fitb <- limma::eBayes(fit)
# })

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
