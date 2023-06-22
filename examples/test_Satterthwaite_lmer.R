library(lmerTest)
library(phylolm)

# From https://bbolker.github.io/mixedmodels-misc/notes/phylog.html

phylo.to.Z <- function(r) {
  ntip <- length(r$tip.label)
  Z <- Matrix(0.0,ncol=length(r$edge.length),nrow=ntip)
  nodes <- (ntip+1):max(r$edge)
  root <- nodes[!(nodes %in% r$edge[,2])]
  for (i in 1:ntip) {
    cn <- i  ## current node
    while (cn != root) {
      ce <- which(r$edge[,2]==cn)   ## find current edge
      Z[i,ce] <- r$edge.length[ce]  ## set Z to branch length
      cn <- r$edge[ce,1]            ## find previous node
    }
  }
  return(Z)
}

#' split a square (block) matrix into component blocks
#' @param M square matrix
#' @param ind indices (0,n1,n2,...) giving the endpoint of each block
split_blkMat <- function(M,ind) {
  res <- list()
  if (length(ind)==1) return(list(M))
  for (i in 1:(length(ind)-1)) {
    v <- (ind[i]+1):ind[i+1]
    res[[i]] <- M[v,v]
  }
  return(res)
}

#' modify reTrms object
#' @param rt a reTrms object
#' @param phylo a phylo object (phylogenetic tree)
#' @param phylonm name of phylogenetic term in model
#' @param phyloZ Z matrix built on branch length
modify_phylo_retrms <- function(rt,phylo,phylonm="phylo",
                                phyloZ=phylo.to.Z(phylo)) {
  ## FIXME: better way to specify phylonm
  ## need to replace Zt, Lind, Gp, flist, Ztlist
  ## we have the same number of parameters (theta, lower),
  ##  same number of obs
  n.edge <- nrow(phylo$edge)
  phylo.pos <- which(names(rt$cnms)==phylonm)
  inds <- c(0,cumsum(sapply(rt$Ztlist,nrow)))
  ## Zt: substitute phylo Z for previous dummy (scalar-intercept) Z
  rt[["Ztlist"]][[phylo.pos]] <- t(phyloZ)
  ## reconstitute Zt from new Ztlist
  rt[["Zt"]] <- do.call(rbind,rt[["Ztlist"]])
  ## Gp: substitute new # random effects (n.edge) for old # (n.phylo)
  Gpdiff <- diff(rt$Gp)  ## old numbers
  Gpdiff_new <- Gpdiff
  Gpdiff_new[phylo.pos] <- n.edge  ## replace
  rt[["Gp"]] <- as.integer(c(0,cumsum(Gpdiff_new)))          ## reconstitute
  ## Lind: replace phylo block with the same element, just more values
  Lind_list <- split(rt[["Lind"]],rep(seq_along(Gpdiff),Gpdiff))
  Lind_list[[phylo.pos]] <- rep(Lind_list[[phylo.pos]][1],n.edge)
  rt[["Lind"]] <- unlist(Lind_list)
  ## Lambdat: replace block-diagonal element in Lambdat with a
  ##   larger diagonal matrix
  Lambdat_list <- split_blkMat(rt[["Lambdat"]],inds)
  Lambdat_list[[phylo.pos]] <- Diagonal(n.edge,1.0)
  rt[["Lambdat"]] <- as(Matrix::.bdiag(Lambdat_list), "dgTMatrix")
  ## flist:
  rt[["flist"]] <- as.list(rt[["flist"]])
  rt[["flist"]][[phylonm]] <- factor(paste0("edge_",seq(n.edge)))
  return(rt)
}

#'
phylo_lmm <- function(formula, data, family, phylo, phyloZ, REML) {
  glmod <- lFormula(formula = formula, data = data, REML = REML)
  glmod$reTrms <- modify_phylo_retrms(glmod$reTrms,phylo,
                                      phylonm="phylo",phyloZ)
  devfun <- do.call(mkLmerDevfun, glmod)
  opt <- optimizeLmer(devfun)
  res <- list(fit = mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr), # mc = match.call(lmer, call("lmer", formula = formula, data = data, REML = REML)))
              devfun = devfun)
  return(res)
}

## Tree
phylo <- read.tree(text = "(((opossum.2:1.359619307e-14,opossum.1:0):1.359619307e-14,opossum.0:0):1,(((cow.2:1.359619307e-14,cow.1:0):1.359619307e-14,cow.0:0):0.6539768865,(((rat.2:1.359619307e-14,rat.1:0):1.359619307e-14,rat.0:0):0.6247450714,(human.1:1.359619307e-14,human.0:0):0.6247450714):0.02923181509):0.3460231135);")
# phylo <- read.tree(text = "(opossum:1,(cow:0.6539768865,(rat:0.6247450714,human:0.6247450714):0.02923181509):0.3460231135);")
n <- length(phylo$tip.label)
id <- phylo$tip.label
design <- paste0("human.", 0:2)
design <- id %in% design
plot(phylo)
tiplabels(pch = 21, col = as.factor(design + 0), bg = as.factor(design + 0))

## Traits
# set.seed(128912)
sim <- phylolm::rTrait(n = 1, phylo, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 1))
traits <- data.frame(species = as.factor(sub("\\.[0-4]$", "", phylo$tip.label)),
                     id = id,
                     phylo = sub("\\.[0-4]$", "", phylo$tip.label),
                     g1 = sim + rnorm(length(sim), 0, sd = sqrt(0.2)),
                     design = as.factor(design + 0))
rownames(traits) <- id
traits$g1[design] <- traits$g1[design] + 0

## Reg
phyloZ <- phylo.to.Z(phylo)

phylo_lmm_fit <- phylo_lmm(g1 ~ design + (1|phylo),
                           data = traits,
                           family = gaussian,
                           phylo = phylo,
                           phyloZ = phyloZ,
                           REML = FALSE)

## phylolm
fit_phylolm <- phylolm::phylolm(g1 ~ design, traits, phylo, model = "BM", measurement_error = TRUE)

summary(phylo_lmm_fit$fit)
summary(fit_phylolm)

all.equal(coef(summary(phylo_lmm_fit$fit))[, 1], coef(summary(fit_phylolm))[, 1])
all.equal(VarCorr(phylo_lmm_fit$fit)$phylo[1]^2, fit_phylolm$sigma2)
all.equal(attr(VarCorr(phylo_lmm_fit$fit), "sc")^2, fit_phylolm$sigma2_error)

## lmerTest
fit_lmerTest <- lmerTest:::as_lmerModLT(phylo_lmm_fit$fit, phylo_lmm_fit$devfun, tol = 1e-08)
aov <- anova(fit_lmerTest)
aov
fit_lmerTest@vcov_varpar
ddf_satterthwaite_BM(fit_phylolm, phylo)
ddf_satterthwaite_sum(fit_phylolm, phylo)

# aov <- anova(fit_lmerTest, ddf = "Kenward-Roger")
# aov
