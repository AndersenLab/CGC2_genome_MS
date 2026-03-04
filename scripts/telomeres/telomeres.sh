#!/bin/bash

genomes=../../processed_data/genomes
output=../../processed_data/genomes

# Identifying telomere repeat count per kb
for fasta in $genomes/c_briggsae*.fa; do
	base=$(basename $fasta)
	trim=${base#c_briggsae.}
	strain=${trim%%.*}

	seqkit locate  --bed -p "TTAGGC" $fasta > $output/${strain}_telomeres.bed
	seqkit fx2tab -n -l $fasta | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' > $output/${strain}_chrom_sizes.bed

	bedtools makewindows -b $output/${strain}_chrom_sizes.bed -w 1000 > $output/${strain}_1kb_windows.bed
	bedtools intersect -a $output/${strain}_1kb_windows.bed -b $output/${strain}_telomeres.bed -c > $output/${strain}_telomeres_binned_1kb.bed
