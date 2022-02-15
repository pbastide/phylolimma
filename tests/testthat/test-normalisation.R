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
