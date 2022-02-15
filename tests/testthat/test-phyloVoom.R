test_that("Comparision with phyloVoom", {
  skip_if_not_installed("phylocompcodeR")
  set.seed(12891026)

  ## Count
  ngenes <- 5000
  nsamples <- 20
  n_diffexp <- 100

  ## Tree
  tree <- ape::rphylo(nsamples, 0.1, 0)
  tree$edge.length <- tree$edge.length / diag(ape::vcv(tree))[1]

  ## PhylocompcodeR
  id.species <- factor(1:nsamples)
  names(id.species) <- tree$tip.label
  mydata <- phylocompcodeR::generateSyntheticData(dataset = "mydata", n.vars = ngenes,
                                                  samples.per.cond = nsamples/2, n.diffexp = n_diffexp,
                                                  id.species = id.species,
                                                  tree = tree,
                                                  prop.var.tree = 1,
                                                  lengths.relmeans = "auto", lengths.dispersions = "auto")

  countMatrix <- mydata@count.matrix
  lengthMatrix <- mydata@length.matrix
  design <- cbind(rep(1, nsamples), mydata@sample.annotations$id.condition)
  rownames(design) <- tree$tip.label

  ## norm factor
  nf <- edgeR::calcNormFactors(countMatrix)
  names(nf) <- tree$tip.label

  # Normalization
  norm.length <- lengthNormalizeRNASeq(countMatrix, lengthMatrix, normalisationFactor = nf,
                                       lengthNormalization = "TPM", dataTransformation = "log2")

  ##############################################################################
  ## voom limma

  # Voom
  norm.voom <- limma::voom(countMatrix, design = NULL,
                           lib.size = colSums(countMatrix / lengthMatrix) * nf * t(lengthMatrix), plot = TRUE)

  # Limma
  fitlimma <- limma::lmFit(norm.voom, design = design)
  # fitlimma <- limma::lmFit(norm.length, design = design)
  fitbayes <- limma::eBayes(fitlimma)
  pvalues <- fitbayes$p.value[, ncol(fitbayes$p.value)]
  adjpvalues <- p.adjust(pvalues, method = 'BH')

  ## FP
  sum(adjpvalues[-(1:n_diffexp)] <= 0.01) / (length(adjpvalues) - n_diffexp)
  ## TP
  sum(adjpvalues[(1:n_diffexp)] <= 0.01) / n_diffexp

  ##############################################################################
  ## limma trend

  # Limma
  fitlimma <- limma::lmFit(norm.length, design = design)
  fitbayes <- limma::eBayes(fitlimma, trend = TRUE)
  pvalues <- fitbayes$p.value[, ncol(fitbayes$p.value)]
  adjpvalues <- p.adjust(pvalues, method = 'BH')

  ## FP
  sum(adjpvalues[-(1:n_diffexp)] <= 0.01) / (length(adjpvalues) - n_diffexp)
  ## TP
  sum(adjpvalues[(1:n_diffexp)] <= 0.01) / n_diffexp

  ##############################################################################
  ## PhyloVoom

  norm.phyloVoom <- phyloVoom(countMatrix, lengthMatrix, normalisationFactor = nf, phy = tree,
                              lengthNormalization = "TPM", dataTransformation = "log2", plot = TRUE)

  fitphylolimma <- phylolmFit(norm.phyloVoom$E, design = design, phy = tree,
                              weights = norm.phyloVoom$weights, measurement_error = TRUE)
  # fitphylolimma <- phylolmFit(norm.length, design = design, phy = tree, measurement_error = TRUE)
  fitphylobayes <- limma::eBayes(fitphylolimma)
  phylopvalues <- fitphylobayes$p.value[, ncol(fitphylobayes$p.value)]
  phyloadjpvalues <- p.adjust(phylopvalues, method = 'BH')

  ## FP
  sum(phyloadjpvalues[-(1:n_diffexp)] <= 0.01) / (length(phyloadjpvalues) - n_diffexp)
  ## TP
  sum(phyloadjpvalues[(1:n_diffexp)] <= 0.01) / n_diffexp

  ##############################################################################
  ## phylolimma trend

  fitphylolimma <- phylolmFit(norm.length, design = design, phy = tree, measurement_error = TRUE)
  fitphylobayes <- limma::eBayes(fitphylolimma, trend = TRUE)
  phylopvalues <- fitphylobayes$p.value[, ncol(fitphylobayes$p.value)]
  phyloadjpvalues <- p.adjust(phylopvalues, method = 'BH')

  ## FP
  sum(phyloadjpvalues[-(1:n_diffexp)] <= 0.01) / (length(phyloadjpvalues) - n_diffexp)
  ## TP
  sum(phyloadjpvalues[(1:n_diffexp)] <= 0.01) / n_diffexp
})
