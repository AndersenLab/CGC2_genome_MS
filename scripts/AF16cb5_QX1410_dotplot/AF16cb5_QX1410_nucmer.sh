#!/bin/bash

# align with nucmer (will spit out a .delta file)
nucmer --maxgap=500 --prefix=AF16cb --coords ../../processed_data/genomes/c_briggsae.QX1410.nanopore.Feb2020.genome.fa ../../processed_data/genomes/c_briggsae.AF16_CB5.genome.fa

# get coordinate file - filter to contain high-quality alignments that are >1kb 
show-coords -r -l -T AF16cb.delta | awk '$5 > 1000' > ../../processed_data/genome_genome_alignments/AF16cb5_QX1410.transformed.tsv

sed -i '1d' ../../processed_data/genome_genome_alignments/AF16cb5_QX1410.transformed.tsv
