#!/bin/bash

chromV_1kb_bed="../../processed_data/scaffolding/CGC2_1kb_windows_chromV.bed"
hifi_bam="../../processed_data/scaffolding/hifi_CGC2_noGaps_chroms_reoriented.bam"
ont_bam="../../processed_data/ONT_CGC2_noGap_chroms_reoriented.bam"
output="../../processed_data/scaffolding"

mosdepth --threads 12 --by $chromV_1kb_bed $output/hifi_chromV_1kb_cov $hifi_bam
mosdepth --threads 12 --by $chromV_1kb_bed $output/ONT_chromV_1kb_cov $ont_bam
