#! /usr/bin/env Rscript

# import relevant packages
library(subSeq)
library(DESeq2)
library(digest)
library(GenomicRanges)

subsample_eset <- function(eset, contrast, methods = c("edgeR", "DESeq2", "voomLimma"),
proportions = 10 ^ seq(-3, 0, .01),
min_reads = 5, presample = NA, nseed = NULL, ...) {
    # technical replicates
    eset <- updateObject(eset) #for old data structures
    tech_group <- eset@colData$technical_replicate_group
    if (!is.null(tech_group)) {
        print("tech rep")
        print(dim(eset@assays[[1]]))
        id <- which(tech_group == '  ')
        if (!(length(id) == 0)) {
            print(id)
            tech_group = as.vector(tech_group)
            print(tech_group)
            tech_group[id] = paste("group_added", id)
            print(tech_group)
        }
        eset <- DESeq2::collapseReplicates(eset, groupby = tech_group)
        print(dim(eset@assays[[1]]))
    }
    
    countmatrix = eset@assays[[1]]
    class = eset@colData$AtlasAssayGroup
    
    # filter for just that contrast
    contrast.v = strsplit(contrast, "_")[[1]]
    counts.contrast = countmatrix[, class %in% contrast.v]
    class.contrast = factor(class[class %in% contrast.v])
    
    # filter countmatrix
    counts.contrast = counts.contrast[rowSums(counts.contrast) >= min_reads, ]

    print(nseed)
    results = subsample(counts.contrast, proportions, method = methods,
    treatment = class.contrast, replications = 1, seed = nseed)
    results
}

# take in an experiment name and a contrast, run subSeq,
# and save to an output folder
args <- commandArgs(trailingOnly = TRUE)
experiment <- args[1]
contrast <- args[2]
presample = as.numeric(args[3])
path_to_esets = args[4]
outfolder = args[5]
trial=as.numeric(args[6])
nseed=as.numeric(args[7])
load(path_to_esets)
print(experiment)
esets <- esets[[experiment]]

print(nseed)
ss <- subsample_eset(eset = esets,
contrast = contrast,
proportions = 10 ^ seq(-3.0, 0, .01),
min_reads = 5,
presample = presample, nseed=nseed)

save(ss, file= paste0(outfolder, "/", experiment, "_", contrast, "_", trial, ".rda"))
