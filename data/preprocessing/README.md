Downloading the Expression Atlas data
====

The following assumes a Linux environment. In this section, work directly in `./atlasdownload`.

To download the latest Atlas Expression data from [here](https://www.ebi.ac.uk/gxa/download.html) (August 28 2018 release), run the following:

    wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/atlas-latest-data.tar.gz

Unpack it into `atlas-latest-data/`.

To create RNASeq_ids.txt, run

    grep htseq2 atlas-latest-data/*/*analysis-methods.tsv | awk '{split($0, a, "/"); print a[2]}' > RNASeq_ids.txt

To copy over .Rdata files, run

    python cp_RS.py RNASeq_ids.txt ../atlas-RNASeq

To get contrasts.csv, run

    ls atlas-latest-data/*/*.pval.bedGraph | awk -F'[\\/|.]' '{OFS = ","; print $2,$4}' > ../contrasts.csv

Processing the ExpressionAtlas data
====

- Change directory to `./scripts/`.

- Update arguments `path_to_rna_seq` and `path_to_contrasts` in the file `preprocess.R`.

- Run `preprocess`:

```R
source("preprocess.R")
preprocess(outfolder = "../../")
```

This will create `contrastsData.rda`, `samplesData.rda`, `esets.rda` in the `data/` directory.

Running the simulations
====

- Initialize variables in `core.R`:
- Change `path_to_contrasts` and `path_to_esets` to the location of where the
`contrasts.rda` and `esets.rda` files
- Submitting jobs to cluster can be controlled with `cluster_args` variable
(see `/superSeq_run.sh` file for settings we used)

- In `/superSeq_run.sh`, or whatever file created by the user, input the path to where `/superSeq_run.R` is
located so the nodes know where to find the file

- In R:

```R
# Load core functions
source("core.R")

# subsample all of the ExpressionAtlas files and save to subsampling_results
# trial is varied to split replications for studies. We used 1:6 (6 total replications). This was done due to cluster constraints.
subsample_all(outfolder = "../../subsampling_results/", trial = 1) 
summary_all(outfolder = "../../subsampling_summary/", infolder = "../..//subsampling_results/", filt_trial = 1)
```
