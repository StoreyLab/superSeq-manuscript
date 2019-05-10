# load relevant data
library(subSeq)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(digest)

# path to contrasts and eset data from the expressionAtlas data
path_to_contrasts <- "../../contrastsData.rda"
path_to_esets <- "../../esets.rda"
path_to_fits <- "../../fits.rda"
# Cluster command for job submission
cluster_args <- "sbatch -p storey ./superSeq.sh"
cluster_args_sum <- "sbatch -p storey ./superSeq_summary.sh"
cluster_args_boot <- "sbatch -p storey ./superSeq_boot.sh"

subsample_all <- function(outfolder, presample = c(1), trial = 1) {
    load(path_to_contrasts)
    dir.create(outfolder, showWarnings = FALSE)
    combos <- broom::inflate(contrastsData, presample = presample)
    combos <- combos %>% group_by(experiment, contrast, presample) %>% mutate(seed = readBin(digest(c(experiment, contrast, presample, trial), raw=TRUE), "integer"))
    cmds <- paste(cluster_args, combos$experiment, combos$contrast,
    combos$presample, path_to_esets, outfolder, trial, combos$seed)
    for (cmd in cmds) {
        cat(cmd, "\n")
        system(cmd)
    }
}

subsample_missing <- function(outfolder, presample = c(1)) {
  # This function exists due to cluster issues
  load(path_to_contrasts)
  dir.create(outfolder, showWarnings = FALSE)
  # all contrasts
  combos <- broom::inflate(contrastsData, presample = presample, trial = 1:6)
  # gather experiments that finished
  fnames <- data.frame(file = list.files(path = "../../subsampling_results/"))
  fnames <- rbind(fnames, data.frame(file = list.files(path = "../../subsampling_results/")))
  out <- fnames %>%
    separate(file, into = c("experiment", "group1", "group2", "trial"), sep = "_", extra = "drop", remove = TRUE) %>%
    separate(trial, into = "trial", sep = ".rda", extra = "drop") %>%
    mutate(contrasts = paste0(group1,"_",group2), trial = as.numeric(trial))
  # extract unfinished jobs
  combos = combos %>% anti_join(out)
  combos <- combos %>% group_by(experiment, contrast, presample, trial) %>% mutate(seed = readBin(digest(c(experiment, contrast, presample, trial), raw=TRUE), "integer"))
  cmds <- paste(cluster_args, combos$experiment, combos$contrast,
                combos$presample, path_to_esets, outfolder, combos$trial, combos$seed)
  for (cmd in cmds) {
    cat(cmd, "\n")
    system(cmd)
  }
}

summary_all <- function(outfolder, infolder, filt_trial = 1) {
    load(path_to_contrasts)
    dir.create(outfolder, showWarnings = FALSE)
    fnames <- data.frame(file = list.files(path = infolder))
    out <- fnames %>%
    separate(file, into = c("experiment", "group1", "group2", "trial"), sep = "_", extra = "drop", remove = TRUE) %>%
    separate(trial, into = "trial", sep = ".rda", extra = "drop") %>%
    mutate(contrasts = paste0(group1,"_",group2), trial = as.numeric(trial), presample = 1)
    out <- out %>% filter(trial == filt_trial) # separate cluser jobs
    cmds <- paste(cluster_args_sum, out$experiment,
    out$contrast, out$presample, outfolder, infolder, out$trial)
    # run subseq calls subsample_contrast
    print(dim(cmds))
    for (cmd in cmds) {
        cat(cmd, "\n")
        system(cmd)
    }
}

import_files <- function(x, infile) {
  print(x)
  read.table(paste0(infile,x))
}

rbind_all2 <- function(infile) {
    all_files <- list.files(infile,recursive = TRUE)
    output <- lapply(all_files, FUN = function(x) import_files(x, infile))
    print("yes")
    out <- dplyr::bind_rows(output)
    saveRDS(out, file = "../../expressionAtlas.rda")
}
