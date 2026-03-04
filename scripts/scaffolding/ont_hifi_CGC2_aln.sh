#!/bin/bash

genome="../../processed_data/genomes/PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed_scaffolds_final.fa"
genome_final="../../processed_data/genomes/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa"
ont="../../data/PB420_ONT_pass.20251006.fastq.gz"
output="../../processed_data/scaffolding"
hifi="../../data/XEAND_20250610_R84050_PL22490-001_1-1-B01_bc2127-bc2127.hifi_reads.bam"


# On scaffolding assembly
minimap2 -ax map-ont -t 24 $genome $ont | samtools sort -@ 8 -o $output/ONT_CGC2.bam
samtools index $output/ONT_CGC2.bam

samtools fastq -@ 24 $hifi \
          | minimap2 -t 24 -ax map-hifi $genome - \
	  | samtools sort -@ 24 -o $output/hifi_to_CGC2.bam

        samtools index -@ 8 $output/hifi_to_CGC2.bam

# On the final assembly after telomric contig re-ordering, gap filling, and chromosome re-orientating 
minimap2 -ax map-ont -t 24 $genome_final $ont | samtools sort -@ 8 -o $output/ONT_CGC2_noGap_chroms_reoriented.bam
samtools index $output/ONT_CGC2_noGap_chroms_reoriented.bam

samtools fastq -@ 24 $hifi \
          | minimap2 -t 24 -ax map-hifi $genome_final - \
          | samtools sort -@ 24 -o $output/hifi_CGC2_noGaps_chroms_reoriented.bam

        samtools index -@ 8 $output/hifi_CGC2_noGaps_chroms_reoriented.bam
