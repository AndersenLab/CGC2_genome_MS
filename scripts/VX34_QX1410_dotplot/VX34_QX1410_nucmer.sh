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
nucmer --maxgap=500 --prefix=VX34_QX1410.contigs --coords ../../processed_data/genomes/c_briggsae.QX1410.nanopore.Feb2020.genome.fa ../../processed_data/genomes/c_briggsae.VX34.nanopore.Feb2020.genome.fa

# get coordinate file - filter to contain high-quality alignments that are >1kb
show-coords -r -l -T VX34_QX1410.contigs.delta | awk -v OFS='\t' '$5 > 1000' > ../../processed_data/genome_genome_alignments/VX34_QX1410.transformed.tsv

sed -i '1d' ../../processed_data/genome_genome_alignments/VX34_QX1410.transformed.tsv
