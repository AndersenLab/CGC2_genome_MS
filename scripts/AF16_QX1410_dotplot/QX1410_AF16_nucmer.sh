#!/bin/bash
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --job-name="mummer"

#activate environment
source activate mummer

# align with nucmer (will spit out a .delta file)
nucmer --maxgap=500 --prefix=AF16_contigs --coords QX1410.ONT.fa AF16.WS276.fa 

# get coordinate file - filter to contain high-quality alignments that are >1kb 
show-coords -r -l -T AF16_contigs.delta | awk '$5 > 1000' > QX1410_AF16_transformed.tsv
