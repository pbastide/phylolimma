## code to prepare `DATASET` dataset goes here

################################################################################
## Read and format count and length data
## From https://github.com/i2bc/InterspeciesDE
## See Bastide et al. 2023
################################################################################
# Read the raw counts
raw_counts <- read.table(file = "https://raw.githubusercontent.com/i2bc/InterspeciesDE/refs/heads/main/data/rawCounts.txt")
raw_counts <- round(raw_counts) # expected RSEM values -> counts
# Read length information
leng <- read.table(file = "https://raw.githubusercontent.com/i2bc/InterspeciesDE/refs/heads/main/data/rawLengths.txt")
# remove NA
raw_counts_noNA <- raw_counts[complete.cases(raw_counts),]
length_noNA <- leng[complete.cases(raw_counts),colnames(raw_counts_noNA)]

################################################################################
## Read and format condition data
################################################################################
condName <- "sights"
# Species condition (1 = blind, 2 = sighted)
states <- read.table(file = "https://raw.githubusercontent.com/i2bc/InterspeciesDE/refs/heads/main/data/states.csv", sep = ",")
condstp <- states$V2
names(condstp) <- states$V1
# Sample Condition
idSpe <- sapply(1:dim(raw_counts)[2],function(i) strsplit(colnames(raw_counts)[i],"_")[[1]][1])
sights <- condstp[idSpe]
names(sights) <- colnames(raw_counts)
sights <- sights - 1
# Format
sights <- factor(sights)
colData <- data.frame(sights)
colData$sights <- factor(colData$sights)
colData$species <- sub("_.*", "", rownames(colData))

################################################################################
## Read and format Tree
################################################################################
library(ape)
## Get the tree
tree <- read.tree(file = "https://raw.githubusercontent.com/i2bc/InterspeciesDE/refs/heads/main/data/crayfish.nodelabels.tre")
tree$node.label <- NULL
# plot(tree)

## Species names
# Match
tree_data_cor <- match(tree$tip.label, colData$species)
data_tree_cor <- match(colData$species, tree$tip.label)
# Species in the tree NOT in data
tree$tip.label[is.na(tree_data_cor)]
# Species in data NOT in the tree
colData$species[is.na(data_tree_cor)]

## Format Tree
# Get rid of species not in data
tree <- drop.tip(tree, tip = tree$tip.label[is.na(tree_data_cor)])
# plot(tree)

################################################################################
## Dataset
################################################################################
crayfish <- list(tree = tree,
                 counts = as.matrix(raw_counts_noNA),
                 lengths = as.matrix(length_noNA),
                 sights = colData)

usethis::use_data(crayfish, overwrite = TRUE)
