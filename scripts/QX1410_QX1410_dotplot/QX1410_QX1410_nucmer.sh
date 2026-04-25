#!/bin/bash

# align with nucmer (will spit out a .delta file)
nucmer --maxgap=500 --prefix=QX_QX --coords ../../processed_data/genomes/c_briggsae.QX1410.nanopore.Feb2020.genome.fa ../../processed_data/genomes/c_briggsae.QX1410.nanopore.Feb2020.genome.fa

# get coordinate file - filter to contain high-quality alignments that are >1kb 
show-coords -r -l -T QX_QX.delta | awk '$5 > 1000' > ../../processed_data/genome_genome_alignments/QX_QX.transformed.tsv

sed -i '1d' ../../processed_data/genome_genome_alignments/QX_QX.transformed.tsv
