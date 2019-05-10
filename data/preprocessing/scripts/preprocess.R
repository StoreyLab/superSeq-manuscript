library(GenomicRanges)
library(dplyr)
library(SummarizedExperiment)
library(plyr)

path_to_rna_seq = "../atlas-RNASeq/"
path_to_contrasts = "../../contrasts.csv"


preprocess <- function(outfolder = "data") {
    library(dplyr)
    # create esets list
    infolder <- path_to_rna_seq
    outfile <- function(f) paste0(outfolder, "/", f, ".rda")
    esets <- list()
    # Combine all summarized experiments into one big list
    for (infile in list.files(infolder, "*.Rdata", full.names = TRUE)) {
        expID <- gsub(".*\\/(E-.*)-atlas.*", "\\1", infile)
        load(infile)
        tmp_eset = experimentSummary[[1]]
        # check for missing count data
        if (ncol(tmp_eset@assays[[1]]) != 0) {
            esets[[expID]] <- tmp_eset
        } else {
            print(infile) # there was only one file with missing count data
        }
    }
    save(esets, file = outfile("esets"))
    
    # Grab various info on experiments
    info_experiment <- function(e) {
        counts <- e@assays[[1]]
        counts.filt <- counts[rowSums(counts) >= 5,]
        d <- data.frame(sample = rownames(e@colData), group = e@colData$AtlasAssayGroup)
        d$depth <- colSums(counts.filt)
        d$genes <- nrow(counts.filt)
        if (!is.null(e@colData$Organism)) {
            d$organism <- e@colData$Organism
        } else {
            d$organism <- e@colData$organism
        }
        d
    }
    
    samplesData <- plyr::ldply(esets, info_experiment, .id = "experiment")
    
    # Load contrasts data and filter of rna-seq studies
    contrastsData <- read.csv(path_to_contrasts,
    header = FALSE, col.names = c("experiment", "contrast"),
    stringsAsFactors = FALSE)
    
    contrastsData <- contrastsData %>%
    filter(experiment %in% names(esets))
    
    extract_filtered <- function(experiment, group1, group2) {
        e <- esets[[experiment[1]]]
        counts <- e@assays[[1]]
        ingroup <- counts[, e@colData$AtlasAssayGroup %in% c(group1, group2)]
        tmp <- ingroup[rowSums(ingroup) >= 5, ]
        data.frame(genes = nrow(tmp), depth = sum(as.numeric(tmp)))
    }
    
    # Filter low count genes
    contrastsData <- contrastsData %>%
    dplyr::select(experiment, contrast) %>%
    tidyr::separate(contrast, c("group1", "group2"), remove = FALSE) %>%
    group_by(experiment, contrast, group1, group2) %>%
    do(extract_filtered(.$experiment, .$group1, .$group2)) %>% 
    ungroup %>%
    mutate(organism = samplesData$organism[match(experiment,
    samplesData$experiment)])
    
    # add replicate-per-group info
    num_replicates <- samplesData %>% dplyr::count(experiment, group)
    
    contrastsData <- contrastsData %>%
    inner_join(num_replicates %>% dplyr::rename(replicates1 = n),
    by = c("experiment", group1 = "group")) %>%
    inner_join(num_replicates %>% dplyr::rename(replicates2 = n),
    by = c("experiment", group2 = "group")) %>%
    tidyr::unite(replicates, replicates1, replicates2,
    sep = "/", remove = FALSE)
    
    save(samplesData, file = outfile("samplesData"))
    save(contrastsData, file = outfile("contrastsData"))
}
