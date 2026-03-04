#!/bin/bash

dir=../../processed_data/scaffolding
nogap_genome="../../processed_data/genomes/CGC2_gapClosing.scaff_seqs"
out=../../processed_data/genomes

# Removing the duplicate haplotig scaffold 7
seqkit grep -v -p "scaffold_7" $nogap_genome > $dir/CGC2_noGaps_noscaff7.fa

# Flipping the chromosome orientations to match QX1410
seqkit grep -f $dir/flip.ids $dir/CGC2_noGaps_noscaff7.fa | seqkit seq -r -p > $dir/CGC2_final_chrom_flipped.tmp.fa

seqkit grep -f $dir/dontFlip.ids $dir/CGC2_noGaps_noscaf7.fa  > $dir/CGC2_finalunflipped.tmp.fa

cat $dir/CGC2_finalunflipped.tmp.fa $dir/CGC2_final_chrom_flipped.tmp.fa > $out/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa

samtools faidx $out/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa

rm $dir/CGC2_final_chrom_flipped.tmp.fa $dir/CGC2_finalunflipped.tmp.fa $dir/CGC2_noGaps_noscaff7.fa
