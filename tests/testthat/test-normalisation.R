test_that("Errors with matrices", {
  set.seed(12891026)

  expect_error(lengthNormalizeRNASeq(1.0, lengthNormalization = "none"),
               "'countMatrix' should be a matrix.")

  ## Counts
  ngenes <- 50
  nsamples <- 10
  countMatrix <- matrix(rbinom(ngenes * nsamples, 10, 0.1), nrow = ngenes)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthNormalization = "none"),
               "Column of count matrix should be named.")
  colnames(countMatrix) <- paste0("s", 1:nsamples)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthNormalization = "none"),
               "Rows of count matrix should be named.")
  rownames(countMatrix) <- paste0("g", 1:ngenes)

  ## Lengths
  expect_error(lengthNormalizeRNASeq(countMatrix),
               "'lengthMatrix' should be a matrix")
  lengthMatrix <- matrix(rbinom(ngenes * nsamples, 1000, 0.5), nrow = nsamples)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix),
               "Count and length matrices should have the same dimension.")
  lengthMatrix <- matrix(rbinom(ngenes * nsamples, 1000, 0.5), nrow = ngenes)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix),
               "Column of length matrix should be named.")
  colnames(lengthMatrix) <- paste0("t", 1:nsamples)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix),
               "Rows of length matrix should be named.")
  rownames(lengthMatrix) <- paste0("q", 1:ngenes)

  ## Names
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix),
               "Count and length matrices should have the same row names.")
  rownames(lengthMatrix) <- paste0("g", 1:ngenes)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix),
               "Count and length matrices should have the same column names.")
  colnames(lengthMatrix) <- paste0("s", 1:nsamples)

  ## Factor
  normalisationFactor <- matrix(1:12, 6)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = normalisationFactor),
               "'normalisationFactor' must be a vector.")
  normalisationFactor <- c(1, 2)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = normalisationFactor),
               "'normalisationFactor' is a vector. Its length should be equal to the number of columns in 'countMatrix'.")
  normalisationFactor <- 1:nsamples
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = normalisationFactor),
               "'normalisationFactor' is a vector. Its names should match the names of columns in 'countMatrix'.")
  names(normalisationFactor) <- paste0("t", 1:nsamples)
  expect_error(lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = normalisationFactor),
               "'normalisationFactor' is a vector. Its names should match the names of columns in 'countMatrix'.")
  names(normalisationFactor) <- paste0("s", 1:nsamples)

})

test_that("Comparision with voom", {
  set.seed(12891026)

  ## Count
  ngenes <- 50
  nsamples <- 10

  countMatrix <- matrix(rbinom(ngenes * nsamples, 10, 0.1), nrow = ngenes)
  colnames(countMatrix) <- paste0("s", 1:nsamples)
  rownames(countMatrix) <- paste0("g", 1:ngenes)

  ## Lengths
  lengthMatrix <- matrix(rbinom(ngenes * nsamples, 1000, 0.5), nrow = ngenes)
  rownames(lengthMatrix) <- paste0("g", 1:ngenes)
  colnames(lengthMatrix) <- paste0("s", 1:nsamples)

  ## norm factor
  nf <- runif(nsamples, 0.0, 1.1)
  names(nf) <- paste0("s", 1:nsamples)

  ## none - log2 vs voom
  norm.length <- lengthNormalizeRNASeq(countMatrix, normalisationFactor = nf,
                                       lengthNormalization = "none", dataTransformation = "log2")
  norm.voom <- limma::voom(countMatrix, design = NULL,
                           lib.size = colSums(countMatrix) * nf)
  expect_equal(norm.length, norm.voom$E)

  ## TPM - log2 vs voom
  norm.length <- lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = nf,
                                       lengthNormalization = "TPM", dataTransformation = "log2")
  norm.voom <- limma::voom(countMatrix, design = NULL,
                           lib.size = colSums(countMatrix / lengthMatrix) * nf * t(lengthMatrix))
  expect_equal(apply(norm.length, 2, sd), apply(norm.voom$E, 2, sd), tolerance = 1e-3)

  ## RPKM - log2 vs voom
  norm.length <- lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = nf,
                                       lengthNormalization = "RPKM", dataTransformation = "log2")
  norm.voom <- limma::voom(countMatrix, design = NULL,
                           lib.size = colSums(countMatrix) * nf * t(lengthMatrix))
  expect_equal(apply(norm.length, 2, sd), apply(norm.voom$E, 2, sd), tolerance = 1e-6)

})


test_that("Comparision with phyloVoom", {
  skip_if_not_installed("phylocompcodeR")
  set.seed(12891026)

  ## Count
  ngenes <- 100
  nsamples <- 20

  ## Tree
  tree <- ape::rphylo(nsamples, 0.1, 0)
  tree$edge.length <- tree$edge.length / diag(ape::vcv(tree))[1]

  # countMatrix <- 10000 + round(t(phylolm::rTrait(ngenes, tree, model = "lambda", parameters = list(sigma2 = 1000000, lambda = 0.8))))
  # colnames(countMatrix) <- tree$tip.label
  # rownames(countMatrix) <- paste0("g", 1:ngenes)
  #
  # ## Lengths
  # lengthMatrix <- round(t(phylolm::rTrait(ngenes, tree, model = "lambda", parameters = list(lambda = 0.5, ancestral.state = 10000))))
  # rownames(lengthMatrix) <- paste0("g", 1:ngenes)
  # colnames(lengthMatrix) <- tree$tip.label

  ## PhylocompcodeR
  id.species <- factor(1:nsamples)
  names(id.species) <- tree$tip.label
  mydata <- phylocompcodeR::generateSyntheticData(dataset = "mydata", n.vars = ngenes,
                                                  samples.per.cond = nsamples/2, n.diffexp = ngenes / 10,
                                                  id.species = id.species,
                                                  tree = tree,
                                                  prop.var.tree = 0.5,
                                                  lengths.relmeans = "auto", lengths.dispersions = "auto")

  countMatrix <- mydata@count.matrix
  lengthMatrix <- mydata@length.matrix
  design <- cbind(rep(1, nsamples), mydata@sample.annotations$id.condition)
  rownames(design) <- tree$tip.label

  ## norm factor
  nf <- edgeR::calcNormFactors(countMatrix)
  names(nf) <- tree$tip.label

  ##############################################################################
  ## none - log2 vs voom

  # Normalization
  norm.length <- lengthNormalizeRNASeq(countMatrix, normalisationFactor = nf,
                                       lengthNormalization = "none", dataTransformation = "log2")

  # Voom
  norm.voom <- limma::voom(countMatrix, design = design,
                           lib.size = colSums(countMatrix) * nf, plot = TRUE)

  expect_equal(norm.length, norm.voom$E)

  # Voom stabilizes weights
  norm.voom2 <- limma::voom(countMatrix, design = design,
                           lib.size = colSums(countMatrix) * nf, plot = TRUE,
                           weights = norm.voom$weights)

  expect_equivalent(norm.voom2$weights[1, ] / mean(norm.voom2$weights[1, ]),
                    rep(1.0, nsamples),
                    tolerance = 0.05)

  # PhyloVoom
  norm.phyloVoom <- phyloVoom(countMatrix, normalisationFactor = nf, phy = tree,
                              design = design,
                              lengthNormalization = "none", dataTransformation = "log2", plot = TRUE)

  expect_equal(norm.phyloVoom$E, norm.voom$E)

  expect_equivalent(norm.voom$weights[1, ] / sum(norm.voom$weights[1, ]),
                    norm.phyloVoom$weights[1, ] / sum( norm.phyloVoom$weights[1, ]),
                    tolerance = 0.01)

  # PhyloVoom stabilizes weights
  norm.phyloVoom2 <- phyloVoom(countMatrix, normalisationFactor = nf, phy = tree,
                               design = design,
                               lengthNormalization = "none", dataTransformation = "log2", plot = TRUE,
                               weights = norm.phyloVoom$weights)

  expect_equivalent(norm.phyloVoom2$weights[1, ] / mean(norm.phyloVoom2$weights[1, ]),
                    rep(1.0, nsamples),
                    tolerance = 0.005)


  ##############################################################################
  ## TPM - log2 vs voom

  # Normalization
  norm.length <- lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = nf,
                                       lengthNormalization = "TPM", dataTransformation = "log2")

  # Voom
  norm.voom <- limma::voom(countMatrix, design = NULL,
                           lib.size = colSums(countMatrix / lengthMatrix) * nf * t(lengthMatrix), plot = TRUE)
  expect_equal(apply(norm.length, 2, sd), apply(norm.voom$E, 2, sd), tolerance = 1e-3)

  # Voom stabilizes weights
  norm.voom2 <- limma::voom(countMatrix, design = design,
                            lib.size = colSums(countMatrix / lengthMatrix) * nf * t(lengthMatrix), plot = TRUE,
                            weights = norm.voom$weights)

  expect_equivalent(norm.voom2$weights[1, ] / mean(norm.voom2$weights[1, ]),
                    rep(1.0, nsamples),
                    tolerance = 0.06)

  # PhyloVoom
  norm.phyloVoom <- phyloVoom(countMatrix, lengthMatrix, normalisationFactor = nf, phy = tree,
                              lengthNormalization = "TPM", dataTransformation = "log2", plot = TRUE)

  expect_equal(norm.phyloVoom$E, norm.length)

  expect_equivalent(norm.voom$weights[1, ] / sum(norm.voom$weights[1, ]),
                    norm.phyloVoom$weights[1, ] / sum( norm.phyloVoom$weights[1, ]),
                    tolerance = 0.05)

  # PhyloVoom stabilizes weights
  norm.phyloVoom2 <- phyloVoom(countMatrix, lengthMatrix, normalisationFactor = nf, phy = tree,
                               design = design,
                               lengthNormalization = "TPM", dataTransformation = "log2", plot = TRUE,
                               weights = norm.phyloVoom$weights)

  expect_equivalent(norm.phyloVoom2$weights[1, ] / mean(norm.phyloVoom2$weights[1, ]),
                    rep(1.0, nsamples),
                    tolerance = 0.04)

  ##############################################################################
  ## RPKM - log2 vs voom

  # Normalization
  norm.length <- lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = nf,
                                       lengthNormalization = "RPKM", dataTransformation = "log2")

  # Voom
  norm.voom <- limma::voom(countMatrix, design = NULL,
                           lib.size = colSums(countMatrix) * nf * t(lengthMatrix))
  expect_equal(apply(norm.length, 2, sd), apply(norm.voom$E, 2, sd), tolerance = 1e-3)

  # Voom stabilizes weights
  norm.voom2 <- limma::voom(countMatrix, design = design,
                            lib.size = colSums(countMatrix) * nf * t(lengthMatrix), plot = TRUE,
                            weights = norm.voom$weights)

  expect_equivalent(norm.voom2$weights[1, ] / mean(norm.voom2$weights[1, ]),
                    rep(1.0, nsamples),
                    tolerance = 0.06)

  # PhyloVoom
  norm.phyloVoom <- phyloVoom(countMatrix, lengthMatrix, normalisationFactor = nf, phy = tree,
                              lengthNormalization = "RPKM", dataTransformation = "log2", plot = TRUE)

  expect_equal(norm.phyloVoom$E, norm.length)

  expect_equivalent(norm.voom$weights[1, ] / sum(norm.voom$weights[1, ]),
                    norm.phyloVoom$weights[1, ] / sum( norm.phyloVoom$weights[1, ]),
                    tolerance = 0.05)

  # PhyloVoom stabilizes weights
  norm.phyloVoom2 <- phyloVoom(countMatrix, lengthMatrix, normalisationFactor = nf, phy = tree,
                               design = design,
                               lengthNormalization = "RPKM", dataTransformation = "log2", plot = TRUE,
                               weights = norm.phyloVoom$weights)

  expect_equivalent(norm.phyloVoom2$weights[1, ] / mean(norm.phyloVoom2$weights[1, ]),
                    rep(1.0, nsamples),
                    tolerance = 0.03)
})
