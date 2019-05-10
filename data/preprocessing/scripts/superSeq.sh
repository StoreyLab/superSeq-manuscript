#!/bin/bash
#SBATCH -J superSeq # a single job name for the array
#SBATCH --qos=1day # which queue
#SBATCH --mem=65000 # mem needed in MB
#SBATCH -n 1 # number of cores
#SBATCH -N 1 # all cores on one machine
#SBATCH --mail-type=END
~/tmp/R-3.4.1/builddir/bin/R -f ~/ajbass/superSeq/simulations/scripts/superSeq_run.R --args $@
