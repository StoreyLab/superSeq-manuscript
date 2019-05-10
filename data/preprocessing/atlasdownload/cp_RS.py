import sys
import os

# usage:
# python cp_RS.py atlas-latest-data/ RNASeq_RData

# intended to be run on a UNIX system (uses grep and cp)
# needs a file RNASeq_ids.txt with the IDs of each RNA-Seq ID in the Atlas
# data. In UNIX, that can be produced with a line like:

# grep htseq2 atlas-latest-data/*/*analysis-methods.tsv | awk '{split($0, a, "/"); print a[2]}' > RNASeq_ids.txt

[input, output] = sys.argv[1:]

if not os.path.exists(output):
    os.makedirs(output)

for l in open("RNASeq_ids.txt"):
    ID = l[:-1]
    
    f = "%s-atlasExperimentSummary.Rdata" % (ID)
    infile = os.path.join("atlas-latest-data", ID, f)
    outfile = os.path.join(output, f)

    # some experiments are baseline; they don't have .Rdata files
    if os.path.exists(infile):
        print ID
        cmd = "cp %s %s" % (infile, outfile)
        os.system(cmd)
