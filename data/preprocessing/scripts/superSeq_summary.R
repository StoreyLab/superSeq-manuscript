#! /usr/bin/env Rscript
library(subSeq)
library(DESeq2)
library(dplyr)
library(qvalue)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
experiment <- args[1]
contrast <- args[2]
presample = as.numeric(args[3])
outfolder = args[4]
infolder = args[5]
trial = as.numeric(args[6])


load(paste0(infolder, "/", experiment, "_", contrast, "_", trial, ".rda"))

ss_out <- summary(ss)

ss_out$experiment <- experiment
ss_out$contrast <- contrast
ss_out$presample <- presample
ss_out$replication <- trial

write.table(ss_out, file = paste0(outfolder, "/", experiment, "_", contrast, "_", trial, ".summary"))