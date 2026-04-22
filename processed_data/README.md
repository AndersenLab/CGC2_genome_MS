# Processed Data

This directory contains the processed datasets used for downstream genomic analyses. Each subdirectory corresponds to a major analysis component. The files are derived from upstream pipelines (see `../scripts/`).

## Directory Structure

processed_data/
- derivative_analysis/
- genome_genome_alignments/
- genomes/
- liftoff/
- marker_liftover/
- orthofinder/
- scaffolding/

## derivative_analysis/

description

### Contains
- `all_derivative.hard.vcf.gz`
- `genotype_matrix.tsv`

## genomes/

description

### Contains
- `AF16_cb5_telomeres_binned_1kb.bed`
- `AF16_telomeres_binned_1kb.bed`
- `c_briggsae.AF16_CB5.genome.fa`
- `c_briggsae.AF16.PRJNA10731.WS276.genome.fa`
- `c_briggsae.AF16.PRJNA10731.WS276.genome.fa.fai`
- `c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa`
- `c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.fai`
- `c_briggsae.QX1410.nanopore.Feb2020.genome.fa`
- `c_briggsae.QX1410.nanopore.Feb2020.genome.fa.fai`
- `c_briggsae.VX34.nanopore.Feb2020.genome.fa`
- `c_briggsae.VX34.nanopore.Feb2020.genome.fa.fai`
- `CGC2_gapClosing.scaff_seqs`
- `CGC2_telomeres_binned_1kb.bed`
- `CGC2_withGaps_telomeres_binned_1kb.bed`
- `PB420.20251025.inbred.withONT.blobFiltered.fa`
- `PB420.20251025.inbred.withONT.blobFiltered.fa.fai`
- `PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed.bin`
- `PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed_scaffolds_final.agp`
- `PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed_scaffolds_final.fa`
- `PB420_telomeric_ctg_reordered_chromIIgapRemoved.fa`
- `QX1410_telomeres_binned_1kb.bed`
- `scaffold_sizes.tsv`
- `VX34_telomeres_binned_1kb.bed`

## liftoff/

description

### Contains
- `AF16toCGC2/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.liftoff.gff`
- `AF16toCGC2/unmapped_features.txt`
- `CGC2toAF16/c_briggsae.AF16.liftoff.gff`
- `CGC2toAF16/unmapped_features.txt`

## marker_liftover/

description

### Contains
- `blastn_results/AF16_marker_hits.tsv.gz`
- `blastn_results/CGC2_marker_hits.tsv.gz`
- `libraries/AF16/c_briggsae.PRJNA10731.WS280.genomic.fa.*`
- `libraries/CGC2/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.*`
- `liftover_results/markerLift_positions_multi_mapping_summ.tsv`
- `liftover_results/markerLift_positions_near_perfect.tsv`
- `liftover_results/markerLift_positions_perfect.tsv`
- `markers/AFHK_indels.clean.tsv`
- `markers/AFHK_indels.tsv`
- `markers/AFHK.primers.fasta`

## orthofinder/

description

### Contains
- `gff/c_briggsae.PRJNA10731.WS280.protein_coding.longest.gff.gz`
- `gff/c_briggsae.QX1410_20250929.csq.longest.gff.gz`
- `gff/CGC2.softMasked.braker.longest.gff.gz`
- `gff/N2.WBonly.WS283.PConly.longest.gff.gz`
- `orthofinder_results/Orthogroups/Orthogroups.GeneCount.tsv`
- `orthofinder_results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt`
- `orthofinder_results/Orthogroups/Orthogroups.tsv`
- `orthofinder_results/Orthogroups/Orthogroups.txt`
- `orthofinder_results/Orthogroups/Orthogroups_UnassignedGenes.tsv`
- `protein/c_briggsae.PRJNA10731.WS280.protein_coding.longest.prot.fa`
- `protein/c_briggsae.QX1410_20250929.csq.longest.protein.fa`
- `protein/CGC2.softMasked.braker.longest.protein.fa`
- `protein/N2.WBonly.WS283.PConly.prot.fa`

## scaffolding/

description

### Contains
- `alignments_sorted.txt`
- `CGC2_1kb_windows_chromV.bed`
- `CGC2_gaps.bed`
- `dontFlip.ids`
- `flip.ids`
- `HiC_contig_coverage.tsv`
- `HiC_PB420_alignment.markedDUP.bam`
- `HiC_PB420_alignment.sorted.bam`
- `hifi_CGC2_noGaps_chroms_reoriented.bam`
- `hifi_chromV_1kb_cov.regions.mosdepth.bed`
- `hifi_to_CGC2.bam`
- `hifi_to_CGC2.bam.bai`
- `ONT_CGC2.bam`
- `ONT_CGC2.bam.bai`
- `ONT_CGC2_noGap_chroms_reoriented.bam`
- `ONT_chromV_1kb_cov.regions.mosdepth.bed`
- `out.hic`
- `out_JBAT.assembly`
- `out_JBAT.hic`
- `scaffold_contig_IDs.tsv`
