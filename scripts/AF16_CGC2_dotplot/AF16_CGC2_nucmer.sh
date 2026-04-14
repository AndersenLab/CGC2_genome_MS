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
nucmer --maxgap=500 --prefix=AF16 --coords ../../processed_data/genomes/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa ../../processed_data/genomes/c_briggsae.AF16.PRJNA10731.WS276.genome.fa

# get coordinate file - filter to contain high-quality alignments that are >1kb 
show-coords -r -l -T AF16.delta | awk '$5 > 1000' > ../../processed_data/genome_genome_alignments/AF16_CGC2.transformed.tsv

sed -i '1d' ../../processed_data/genome_genome_alignments/AF16_CGC2.transformed.tsv
