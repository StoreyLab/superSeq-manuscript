The following assumes a Linux environment.

To download the latest Atlas Expression data from [here](https://www.ebi.ac.uk/gxa/download.html), run the following (you might want to check that the link hasn't changed):

    wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/atlas-latest-data.tar.gz

Unpack it into atlas-latest-data.

To create RNASeq_ids.txt, run

    grep htseq2 atlas-latest-data/*/*analysis-methods.tsv | awk '{split($0, a, "/"); print a[2]}' > RNASeq_ids.txt

To copy over .Rdata files, run

    python cp_RS.py RNASeq_ids.txt ../atlas-RNASeq

To get contrasts.csv, run

    ls atlas-latest-data/*/*.pval.bedGraph | awk -F'[\\/|.]' '{OFS = ","; print $2,$4}' > ../contrasts.csv
