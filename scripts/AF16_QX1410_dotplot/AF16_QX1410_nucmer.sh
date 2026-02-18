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
nucmer --maxgap=500 --prefix=AF16 --coords ../../raw_data/genomes/c_briggsae.QX1410_nanopore.Feb2020.genome.fa ../../raw_data/genomes/c_briggsae.PRJNA10731.WS276.AF16.genome.fa

# get coordinate file - filter to contain high-quality alignments that are >1kb 
show-coords -r -l -T AF16.delta | awk '$5 > 1000' > ../../processed_data/genome_genome_alignments/AF16_QX1410.transformed.tsv
