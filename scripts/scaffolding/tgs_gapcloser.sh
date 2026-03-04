#!/bin/bash

asm="../../processed_data/genomes/PB420_telomeric_ctg_reordered_chromIIgapRemoved.fa"
ont="../../data/PB420_ONT_pass.20251006.fastq.gz"
ont_fa="${ont%.fastq.gz}.fa"
sr1="../../data/PB420_POOLRET57_1N_1R.fq.gz"
sr2="../../data/PB420_POOLRET57_1N_2R.fq.gz"
SR="../../data/PB420_SR_merged.fastq.gz"

seqkit fq2fa $ont > $ont_fa
cat $sr1 $sr2 > $SR

tgsgapcloser \
	--scaff $asm \
	--reads $ont_fa \
	--pilon $pilon_jar \
	--ngs $SR \
	--output CGC2_gapClosing \
	--thread 24 \
	--samtools $samtools_install \
	--java $java
