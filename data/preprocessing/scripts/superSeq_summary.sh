#!/bin/bash
#SBATCH -J superSeq # a single job name for the array
#SBATCH --qos=1day # which queue
#SBATCH --mem=50000 # mem needed in MB
#SBATCH -n 1 # number of cores
#SBATCH -N 1 # all cores on one machine
#SBATCH --mail-type=END
~/tmp/R-3.4.1/builddir/bin/R -f /Genomics/storeyscratch/storeylab/ajbass/ajbass/superSeq/simulations/scripts/superSeq_summary.R --args $@
