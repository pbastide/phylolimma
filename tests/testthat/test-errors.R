test_that("Errors with species names", {
  set.seed(12891026)
  ## Tree
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)
  mat_tree <- ape::vcv(tree)

  ## data
  ngenes <- 50
  # wrong dimension
  y_data <- matrix(rnorm(ngenes*ntips), ncol = ngenes)
  expect_error(phylolmFit(y_data, phy = tree),
               "not equal to array extent")
  expect_error(checkParamMatrix(y_data, "data", tree),
               "`data` should have as many columns as the number of taxa in the tree.")
  # no names
  y_data <- matrix(rnorm(ngenes*ntips), nrow = ngenes)
  expect_error(phylolmFit(y_data, phy = tree),
               "`expression matrix` and/or the tips of the phylogeny are not named.")
  # wrong order
  colnames(y_data) <- sample(tree$tip.label, ntips)
  expect_warning(res1 <- phylolmFit(y_data, phy = tree),
                 "expression matrix` was not sorted in the correct order")
  y_data <- y_data[, match(tree$tip.label, colnames(y_data)), drop = FALSE]
  res2 <- phylolmFit(y_data, phy = tree)
  expect_equal(res1, res2)
  # wrong names
  colnames(y_data)[1] <- "moustache"
  expect_error(phylolmFit(y_data, phy = tree),
               "Species 't1' are in the tree but not in expression matrix.\n  Species 'moustache' is in expression matrix but not in the tree.")
  # Correct name and order
  colnames(y_data) <- tree$tip.label
  # wrong names tree
  tree_wrong <- tree
  tree_wrong$tip.label[1] <- "moustache"
  expect_error(phylolmFit(y_data, phy = tree_wrong),
               "Species 'moustache' are in the tree but not in expression matrix.\n  Species 't1' is in expression matrix but not in the tree.")


  ## Design
  design <- matrix(1, nrow = ntips, ncol = 2)
  design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
  # wrong dimension
  expect_error(phylolmFit(y_data, design = t(design), phy = tree),
               "row dimension of design doesn't match column dimension of data object")
  expect_error(checkParamMatrix(t(design), "design", tree, TRUE),
               "`design` should have as many rows as the number of taxa in the tree.")
  # no names
  expect_error(phylolmFit(y_data, design = design, phy = tree),
               "`design matrix` and/or the tips of the phylogeny are not named.")
  # wrong order
  rownames(design) <- sample(tree$tip.label, ntips)
  expect_warning(res1 <- phylolmFit(y_data, design = design, phy = tree),
                 "`design matrix` was not sorted in the correct order")
  design <- design[match(tree$tip.label, rownames(design)), , drop = FALSE]
  res2 <- phylolmFit(y_data, design = design, phy = tree)
  expect_equal(res1, res2)
  # wrong names
  rownames(design)[1] <- "moustache"
  expect_error(phylolmFit(y_data, design = design, phy = tree),
               "Species 't1' are in the tree but not in design matrix.\n  Species 'moustache' are in design matrix but not in the tree.")
})

test_that("Unused parameters", {
  set.seed(12891026)
  ## Tree
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)
  mat_tree <- ape::vcv(tree)
  ## data
  ngenes <- 50
  # wrong dimension
  y_data <- matrix(rnorm(ngenes*ntips), nrow = ngenes)
  colnames(y_data) <- tree$tip.label

  expect_error(phylolmFit(y_data, phy = tree, ndups = 2), "'ndups' can only be '1'")
  expect_error(phylolmFit(y_data, phy = tree, spacing = 2), "'spacing' can only be '1'")
  # expect_error(phylolmFit(y_data, phy = tree, weights = rep(0.5, ntips)), "'weights' can only be 'null'")
  expect_error(phylolmFit(y_data, phy = tree, method = "robust"), "'method' can only be 'ls'")
  expect_error(phylolmFit(y_data, phy = tree, correlation = 0.5), "'correlation' is not used")
  expect_error(phylolmFit(y_data, phy = tree, block = sample(c(0, 1), ntips, replace = TRUE)), "'block' can only be 'null'")

  y <- limma::getEAWP(y_data)
  expect_error(phylolmFit(y, phy = tree), "'object' must be a matrix.")
})
