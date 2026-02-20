#!/bin/bash

output=../telomeres

# For CGC2
#genome="/vast/eande106/data/c_briggsae/genomes/CGC2/Jan2026/c_briggsae.CGC2_HiFi.ONT.HiC.20260122.fa"

#seqkit locate  --bed -p "TTAGGC" $genome > $output/CGC2_telomeres.bed

#seqkit fx2tab -n -l $genome | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' > $output/CGC2_chrom_sizes.bed

#bedtools makewindows -b $output/CGC2_chrom_sizes.bed -w 1000 > $output/CGC2_1kb_windows.bed

#bedtools intersect -a $output/CGC2_1kb_windows.bed -b $output/CGC2_telomeres.bed -c > $output/CGC2_telomeres_binned_1kb.bed



# For VX34
#genome="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/VX34.fa"

#seqkit locate  --bed -p "TTAGGC" $genome > $output/VX_telomeres.bed

#seqkit fx2tab -n -l $genome | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' > $output/VX_chrom_sizes.bed

#bedtools makewindows -b $output/VX_chrom_sizes.bed -w 1000 > $output/VX_1kb_windows.bed

#bedtools intersect -a $output/VX_1kb_windows.bed -b $output/VX_telomeres.bed -c > $output/VX_telomeres_binned_1kb.bed


# For QX1410
#genome="/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa"

#seqkit locate  --bed -p "TTAGGC" $genome > $output/QX_telomeres.bed

#seqkit fx2tab -n -l $genome | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' > $output/QX_chrom_sizes.bed

#bedtools makewindows -b $output/QX_chrom_sizes.bed -w 1000 > $output/QX_1kb_windows.bed

#bedtools intersect -a $output/QX_1kb_windows.bed -b $output/QX_telomeres.bed -c > $output/QX_telomeres_binned_1kb.bed




# For AF16
genome="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/scripts_NOTGITHUB/AF16.WormBase.fa"

#seqkit locate  --bed -p "TTAGGC" $genome > $output/af_telomeres.bed

#seqkit fx2tab -n -l $genome | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' > $output/af_chrom_sizes.bed

bedtools makewindows -b $output/af_chrom_sizes.bed -w 1000 > $output/af_1kb_windows.bed

bedtools intersect -a $output/af_1kb_windows.bed -b $output/af_telomeres.bed -c > $output/af_telomeres_binned_1kb.bed







