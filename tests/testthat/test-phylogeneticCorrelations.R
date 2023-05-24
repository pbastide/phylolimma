test_that("phylogeneticCorrelations - BM", {
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

  #################################################################################################
  ## No measurement error

  measurement_error <- FALSE
  ## Fit Phylo
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = model,
                              measurement_error = measurement_error,
                              use_consensus = FALSE)
  ## Fit Phylo Consensus
  resPhyloLmFitCons <- phylolmFit(y_data, design = design, phy = tree,
                                  model = model,
                                  measurement_error = measurement_error,
                                  use_consensus = TRUE)

  ## The two objects should be equal
  expect_equal(resPhyloLmFit, resPhyloLmFitCons)

  #################################################################################################
  ## With measurement error

  measurement_error <- TRUE
  ## Fit Phylo
  resPhyloLmFit <- phylolmFit(y_data, design = design, phy = tree,
                              model = model,
                              measurement_error = measurement_error,
                              use_consensus = FALSE)
  ## Fit Phylo Consensus
  resPhyloLmFitCons <- phylolmFit(y_data, design = design, phy = tree,
                                  model = model,
                                  measurement_error = measurement_error,
                                  use_consensus = TRUE)

  ## Test names and dimensions
  expect_equal(colnames(resPhyloLmFit$coefficients), colnames(resPhyloLmFitCons$coefficients))
  expect_equal(dim(resPhyloLmFit$coefficients), dim(resPhyloLmFitCons$coefficients))
  expect_equal(length(resPhyloLmFit$df.residual), length(resPhyloLmFitCons$df.residual))
  expect_equal(length(resPhyloLmFit$sigma), length(resPhyloLmFitCons$sigma))
  expect_equal(colnames(resPhyloLmFit$stdev.unscaled), colnames(resPhyloLmFitCons$stdev.unscaled))
  expect_equal(dim(resPhyloLmFit$stdev.unscaled), dim(resPhyloLmFitCons$stdev.unscaled))

  ## ebayes
  fitphy <- eBayes(resPhyloLmFit)
  fitphycons <- eBayes(resPhyloLmFitCons)

  ## Test names and dimensions
  expect_equal(colnames(fitphy$p.value), colnames(fitphycons$p.value))
  expect_equal(dim(fitphy$p.value), dim(fitphycons$p.value))

  #################################################################################################
  ## Pagel lambda

  ## Fit Phylo Consensus
  resPhyloLmFitConsLambda <- phylolmFit(y_data, design = design, phy = tree,
                                        model = "lambda",
                                        measurement_error = FALSE,
                                        use_consensus = TRUE)

  expect_equal(resPhyloLmFitConsLambda@.Data[which(!names(resPhyloLmFitConsLambda) %in% c("modelphy", "measurement_error"))],
               resPhyloLmFitCons@.Data[which(!names(resPhyloLmFitConsLambda) %in% c("modelphy", "measurement_error"))],
               tol = 1e-5)


})

test_that("phylogeneticCorrelations - separate call", {
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

  ## checks same fit
  for(model in c("BM", "lambda", "OUfixedRoot", "delta")) {
    for (measurement_error in c(FALSE, TRUE)) {
      if (model == "lambda" && measurement_error) {
        expect_error(
          phylolmFit(y_data, design = design, phy = tree,
                     model = model,
                     measurement_error = measurement_error,
                     use_consensus = TRUE),
        "the lambda transformation and measurement error cannot be used together: they are not distinguishable")
        expect_error(
          phylogeneticCorrelations(y_data, design = design, phy = tree,
                     model = model,
                     measurement_error = measurement_error,
                     use_consensus = TRUE),
          "the lambda transformation and measurement error cannot be used together: they are not distinguishable")
      } else {
        res1 <- phylolmFit(y_data, design = design, phy = tree,
                           model = model,
                           measurement_error = measurement_error,
                           use_consensus = TRUE,
                           trim = 0.25)
        pc <- phylogeneticCorrelations(y_data, design = design, phy = tree,
                                       model = model,
                                       measurement_error = measurement_error,
                                       trim = 0.25)
        res2 <- phylolmFit(y_data, design = design, phy = tree,
                           model = model,
                           measurement_error = measurement_error,
                           use_consensus = TRUE,
                           consensus_tree = pc)
        expect_equal(res1, res2)
      }
    }
  }

  ## Same fit but different trim gives different result
  res1 <- phylolmFit(y_data, design = design, phy = tree,
                     model = model,
                     measurement_error = measurement_error,
                     use_consensus = TRUE)
  pc <- phylogeneticCorrelations(y_data, design = design, phy = tree, model = model, measurement_error = measurement_error, trim = 0.9)
  res2 <- phylolmFit(y_data, design = design, phy = tree,
                     model = model,
                     measurement_error = measurement_error,
                     use_consensus = TRUE,
                     consensus_tree = pc)
  expect_true(res1$coefficients[1, 1] != res2$coefficients[1, 1])



})

test_that("phylogeneticCorrelations - eBayes", {
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

  ## Fit Phylo Consensus
  fit1 <- phylolmFit(y_data, design = design, phy = tree,
                     model = "OUfixedRoot",
                     measurement_error = TRUE,
                     use_consensus = TRUE)
  fit1eb <- eBayes(fit1)
  fit1ebtrend <- eBayes(fit1, trend = TRUE)

  ## phyCor
  pc <- phylogeneticCorrelations(y_data, design = design, phy = tree,
                                 model = "OUfixedRoot",
                                 measurement_error = TRUE)
  fit2 <- phylolmFit(y_data, design = design, phy = tree,
                     model = "OUfixedRoot",
                     measurement_error = TRUE,
                     use_consensus = TRUE,
                     consensus_tree = pc)
  fit2eb <- eBayes(fit2)
  fit2ebtrend <- eBayes(fit2, trend = TRUE)

  expect_equal(fit1, fit2)
  expect_equal(fit1eb, fit2eb)
  expect_equal(fit1ebtrend, fit2ebtrend)
})


test_that("phylogeneticCorrelations - Errors", {
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

  #################################################################################################
  ## Errors

  # matrix of data
  expect_error(phylogeneticCorrelations(y_data[1, ], design = design, phy = tree,
                                        model = model,
                                        measurement_error = TRUE,
                                        trim = 0.25),
               "must be a matrix.")
  expect_error(phylolmFit(y_data[1, ], design = design, phy = tree,
                          model = model,
                          measurement_error = TRUE,
                          use_consensus = TRUE,
                          trim = 0.25),
               "must be a matrix.")

  # no weights
  expect_error(phylogeneticCorrelations(y_data, design = design, phy = tree,
                                        model = model,
                                        measurement_error = TRUE,
                                        trim = 0.25, weights = 1:ntips / sum(1:ntips)),
               "weights are not allowed with the phylogenetic regression")
  expect_error(phylolmFit(y_data[1, ], design = design, phy = tree,
                          model = model,
                          measurement_error = TRUE,
                          use_consensus = TRUE,
                          trim = 0.25, weights = 1:ntips / sum(1:ntips)),
               "'weights' can only be 'null' in 'phylolmFit'")

  # numeric design
  designError <- design
  designError[, "condition"] <- c("A", "B")[designError[, "condition"] + 1]
  expect_error(phylogeneticCorrelations(y_data, design = designError, phy = tree,
                                        model = model,
                                        measurement_error = TRUE,
                                        trim = 0.25),
               "design must be a numeric matrix")
  expect_error(phylolmFit(y_data, design = designError, phy = tree,
                          model = model,
                          measurement_error = TRUE,
                          use_consensus = TRUE,
                          trim = 0.25),
               "design must be a numeric matrix")

  # identifiable design
  designError <- design
  designError <- cbind(designError, 1 - designError[, "condition"])
  expect_error(phylogeneticCorrelations(y_data, design = designError, phy = tree,
                                        model = model,
                                        measurement_error = TRUE,
                                        trim = 0.25),
               "Coefficients not estimable:")
  expect_error(phylolmFit(y_data, design = designError, phy = tree,
                          model = model,
                          measurement_error = TRUE,
                          use_consensus = TRUE,
                          trim = 0.25),
               "Coefficients not estimable:")

  # tree
  treeError <- tree
  attr(treeError, "class") <- NULL
  expect_error(phylogeneticCorrelations(y_data, design = design, phy = treeError,
                                        model = model,
                                        measurement_error = TRUE,
                                        trim = 0.25),
               "must be of class 'phylo'")
  expect_error(phylolmFit(y_data, design = design, phy = treeError,
                          model = model,
                          measurement_error = TRUE,
                          use_consensus = TRUE,
                          trim = 0.25),
               "must be of class 'phylo'")
})
