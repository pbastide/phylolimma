test_that("phylogeneticCorrelations - star tree", {
  skip_if_not_installed("statmod")
  skip_if_not_installed("phytools")
  set.seed(12891026)
  ## Star tree
  ntips <- 50
  tree <- ape::stree(ntips, "star")
  tree$edge.length <- rep(1.0, ntips)
  ## Replicates
  r <- 2
  ids <- as.vector(sapply(1:r, function(i) paste0(tree$tip.label, "_", i)))
  traits <- data.frame(species = sub("\\_.", "", ids), ids = ids)
  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "ids")

  ## data
  ngenes <- 40
  y_data <- t(phylolm::rTrait(ngenes, tree_rep, model = "lambda", parameters = list(lambda = 0.8)))
  ## Design
  design <- matrix(1, nrow = ntips, ncol = 2)
  design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
  colnames(design) <- c("(Intercept)", "condition")
  design <- design[rep(seq_len(nrow(design)), each = r), ]
  rownames(design) <- tree_rep$tip.label
  y_data[design[, "condition"]] <- y_data[design[, "condition"]] + 2

  ## Phylogenetic Correlations
  phycor <- phylogeneticCorrelations(y_data, design = design, phy = tree_rep,
                                     model = "BM",
                                     measurement_error = TRUE,
                                     trim = 0.15)

  phycor2 <- phylogeneticCorrelations(y_data, design = design, phy = tree_rep,
                                      model = "lambda",
                                      measurement_error = FALSE,
                                      trim = 0.15)

  ## Duplicated Correlations
  ducor <- limma::duplicateCorrelation(y_data, design = design, ndups = 1,
                                       block = sub("_[1-9]", "", tree_rep$tip.label))

  expect_equal(phycor$params$lambda_error, ducor$cor, tolerance = 1e-2)
  expect_equal(phycor2$params$lambda, ducor$cor, tolerance = 1e-2)

  ## Fits
  fitlimma <- limma::lmFit(y_data, design = design,
                           correlation = ducor$consensus,
                           block = sub("_[1-9]", "", tree_rep$tip.label))
  fitphyloDE <- phylolmFit(y_data, design = design,
                           phy = tree_rep, model = 'BM',
                           measurement_error = TRUE,
                           use_consensus = TRUE, consensus_tree = phycor)

  expect_equal(fitlimma$sigma, fitphyloDE$sigma, tolerance = 1e-5)
  expect_equal(fitlimma$coefficients, fitphyloDE$coefficients, tolerance = 1e-5)
  expect_equal(fitlimma$design, fitphyloDE$design)

  ## eBayes
  fitbayes <- limma::eBayes(fitlimma, trend = FALSE)
  fitphylobayes <- limma::eBayes(fitphyloDE, trend = FALSE)

  expect_equal(fitbayes$t, fitphylobayes$t, tolerance = 1e-5)
  expect_equal(fitbayes$p.value, fitphylobayes$p.value, tolerance = 1e-5)

  expect_equal(fitbayes$cov.coefficients, fitphylobayes$cov.coefficients, tolerance = 1e-5)
  expect_equal(fitbayes$F, fitphylobayes$F, tolerance = 1e-5)
  expect_equal(fitbayes$F.p.value, fitphylobayes$F.p.value, tolerance = 1e-5)

  ## treat
  fittreat <- limma::treat(fitlimma, robust = TRUE)
  fitphylotreat <- limma::treat(fitphyloDE, robust = TRUE)

  expect_equal(fittreat$t, fitphylotreat$t, tolerance = 1e-5)
  expect_equal(fittreat$p.value, fitphylotreat$p.value, tolerance = 1e-5)

  ## decideTests
  fitdt <- limma::decideTests(fitlimma, method = "global")
  fitphylodt <- limma::decideTests(fitphyloDE, method = "global")
  expect_equal(fitdt@.Data, fitphylodt@.Data)

  fitdt <- limma::decideTests(fitlimma, method = "nestedF")
  fitphylodt <- limma::decideTests(fitphyloDE, method = "nestedF")
  expect_equal(fitdt@.Data, fitphylodt@.Data)

  fitdt <- limma::decideTests(fitlimma, method = "hierarchical", p.value = 0.5)
  fitphylodt <- limma::decideTests(fitphyloDE, method = "hierarchical", p.value = 0.5)
  expect_equal(fitdt@.Data, fitphylodt@.Data)

})

test_that("phylogeneticCorrelations - Convergence issues", {
  skip_if_not_installed("statmod")
  skip_if_not_installed("phytools")
  set.seed(12891026)
  ## Star tree
  ntips <- 50
  tree <- ape::stree(ntips, "star")
  tree$edge.length <- rep(1.0, ntips)
  ## Replicates
  r <- 2
  ids <- as.vector(sapply(1:r, function(i) paste0(tree$tip.label, "_", i)))
  traits <- data.frame(species = sub("\\_.", "", ids), ids = ids)
  tree_rep <- addReplicatesOnTree(tree, traits, species = "species", id = "ids")

  ## data
  ngenes <- 40
  y_data <- t(phylolm::rTrait(ngenes - 1, tree_rep, model = "lambda", parameters = list(lambda = 0.8)))
  y_data_dege <- phylolm::rTrait(1, tree, model = "lambda", parameters = list(lambda = 1))
  y_data_dege <- rep(y_data_dege, each = r)
  y_data <- rbind(y_data, y_data_dege)
  ## Design
  design <- matrix(1, nrow = ntips, ncol = 2)
  design[sample(1:ntips, floor(ntips / 2)), 2] <- 0
  colnames(design) <- c("(Intercept)", "condition")
  design <- design[rep(seq_len(nrow(design)), each = r), ]
  rownames(design) <- tree_rep$tip.label
  y_data[design[, "condition"]] <- y_data[design[, "condition"]] + 2

  ## Phylogenetic Correlations
  phycor <- phylogeneticCorrelations(y_data, design = design, phy = tree_rep,
                                     model = "lambda",
                                     measurement_error = FALSE,
                                     trim = 0.15)

  phycor2 <- phylogeneticCorrelations(y_data, design = design, phy = tree_rep,
                                      lower.bound = list(sigma2_error = 1e-20),
                                      model = "BM",
                                      measurement_error = TRUE,
                                      trim = 0.15)

  ## Duplicated Correlations
  ducor <- limma::duplicateCorrelation(y_data, design = design, ndups = 1,
                                       block = sub("_[1-9]", "", tree_rep$tip.label))

  expect_equal(phycor$params$lambda, ducor$cor, tolerance = 1e-2)
  expect_equal(phycor2$params$lambda_error, phycor$params$lambda, tolerance = 1e-2)

})
