#!/bin/bash

output=../alignment
asm="../../processed_data/genomes/PB420.20251025.inbred.withONT.blobFiltered.fa"
r1="../../data/pb420_hic_r1.fastq.gz"
r2="../../data/pb420_hic_r2.fastq.gz"

# Alignment of Hi-C seq data to contiguous assembly
$bwa_mem2 index $asm
$bwa_mem2 mem -t 24 -5 -S -P $asm $r1 $r2 | samtools view -@ 24 -b | samtools sort -@ 24 -o $output/HiC_PB420_alignment.sorted.bam 
