#!/bin/bash

workdir=../../processed_data/scaffolding

# Adding RG labels
picard AddOrReplaceReadGroups \
  I=$workdir/HiC_PB420_alignment.sorted.bam \
  O=$workdir/HiC_PB420_alignment.sorted.RG.bam \
  RGID=PB420_HiC \
  RGLB=HiC \
  RGPL=ILLUMINA \
  RGPU=PB420_HiC \
  RGSM=PB420 \
  VALIDATION_STRINGENCY=SILENT

# Marking duplicate Hi-C reads
picard -Xmx150g MarkDuplicates \
  I=$workdir/HiC_PB420_alignment.sorted.RG.bam \
  O=$workdir/HiC_PB420_alignment.markedDUP.bam \
  M=$workdir/HiC_PB420_alignment.Picard_dup_metrics.txt \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true \
  TMP_DIR=$TMPDIR
