#!/bin/bash

bam="../../processed_data/scaffolding//HiC_PB420_alignment.markedDUP.bam"
asm="../../processed_data/genomes/PB420.20251025.inbred.withONT.blobFiltered.fa"
outdir="../../processed_data/genomes"

$yahs $asm $bam -o $outdir/PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffolded
