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

## genome_genome_alignments/

This directory contains TSVs of genome-genome alignment coordinates that are output from nucmer and used to create dotplots. 

### Contains
- `AF16cb5_QX1410.transformed.tsv`
	- Genome-genome alignment coordinates of AF16 cb5 to QX1410.
- `AF16_CGC2.transformed.tsv`
	- Genome-genome alignment coordinates of AF16 cb4 to CGC2.
- `AF16_QX1410.transformed.tsv`
	- Genome-genome alignment coordinates of AF16 cb4 to QX1410.
- `CGC2_QX1410.transformed.tsv`
	- Genome-genome alignment coordinates of CGC2 to QX1410.
- `CGC2_withGaps_QX1410.transformed.tsv`
	- Genome-genome alignment coordinates of CGC2 (before gap closing) to QX1410.
- `QX_QX.transformed.tsv`
	- Genome-genome alignment coordinates of QX1410 to itself.
- `scaffold7_CGC2.transformed.tsv`
	- Genome-genome alignment coordinates of CGC2 scaffold 7 to CGC2
- `VX34_QX1410.transformed.tsv`
	- Genome-genome alignment coordinates of VX34 to QX1410.

## genomes/

This directory contains *C. briggsae* reference genomes, their indices, and quantitative metrics of their assemblies such as where telomeric repeats are found. Additionally, this directory contains the contiguous CGC2 assembly and all intermediate stages of the genome after scaffolding and the manual curation process.

### Contains
- `AF16_cb5_telomeres_binned_1kb.bed`
	- Counts of telomeric repeats per kb in the AF16 cb5 genome.
- `AF16_telomeres_binned_1kb.bed`
	- Counts of telomeric repeats per kb in the AF16 cb4 genome.
- `c_briggsae.AF16_CB5.genome.fa`
	- The AF16 cb5 genome retrieved from NCBI under accession PRJNA917437.
- `c_briggsae.AF16.PRJNA10731.WS276.genome.fa`
	- The AF16 cb4 genome retrieved from NCBI.
- `c_briggsae.AF16.PRJNA10731.WS276.genome.fa.fai`
	- FASTA index of the AF16 cb4 genome.
- `c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa`
	- The final CGC2 genome assembly available on NCBI under accession PRJNA1451919.
- `c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.fai`
	- FASTA index of the final CGC2 genome assembly.
- `c_briggsae.QX1410.nanopore.Feb2020.genome.fa`
	- The QX1410 genome retrieved from NCBI under accession PRJNA784955.
- `c_briggsae.QX1410.nanopore.Feb2020.genome.fa.fai`
	- FASTA index of the QX1410 genome.
- `c_briggsae.VX34.nanopore.Feb2020.genome.fa`
	- The VX34 genome retrieved from NCBI under accession PRJNA784955.
- `c_briggsae.VX34.nanopore.Feb2020.genome.fa.fai`
	- FASTA index of the VX34 genome.
- `CGC2_gapClosing.scaff_seqs` This is a large file not present in this repo.
	- Sequence of CGC2 scaffolds. This file is produced by `../../scripts/scaffolding/tgs_gapcloser.sh`
- `CGC2_telomeres_binned_1kb.bed`
	- Counts of telomeric repeats per kb in the CGC2 genome.
- `CGC2_withGaps_telomeres_binned_1kb.bed`
	- Counts of telomeric repeats per kb in the CGC2 genome before gap closing and re-ordering of contigs on the left-end of chromosome V.
- `PB420.20251025.inbred.withONT.blobFiltered.fa` Large file. Available publicly on AWS.
	- Contiguous assembly of CGC2 before scaffolding.
- `PB420.20251025.inbred.withONT.blobFiltered.fa.fai`
	- FASTA index of contiguous CGC2 assembly.
- `PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffolded.bin` This is a large file not present in this repo.
	- Output file from scaffolding with YaHS that is required for generating a Hi-C contact map with Juicer. This file is produced from `../../scripts/scaffolding/yahs.sh`.
- `PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed_scaffolds_final.agp`
	- Output file from scaffolding with YaHS indicating contig placement into scaffolds and the loci of gaps.
- `PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed_scaffolds_final.fa` This is a large file not present in this repo.
	- Scaffolded CGC2 genome assembly output from YaHS. This file is produced from `../../scripts/scaffolding/yahs.sh`.
- `QX1410_telomeres_binned_1kb.bed`
	- Counts of telomeric repeats per kb in the QX1410 genome.
- `scaffold_sizes.tsv`
	- The sizes of CGC2 scaffolds.
- `VX34_telomeres_binned_1kb.bed`
	- Counts of telomeric repeats per kb in the VX34 genome.

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

This directory contains files that are required for scaffolding with YaHS, Hi-C map creation, and manual curation of CGC2.

### Contains
- `alignments_sorted.txt` This is a large file not present in this repo.
	- Ouput from Juicer of Hi-C contact pairs required for creating a Hi-C contact map. This file is produced from `../../scripts/scaffolding/hiC_contactMap.sh`
- `CGC2_1kb_windows_chromV.bed`
	- CGC2 chromosome V broken up into 1 kb bins.
- `CGC2_gaps.bed`
	- Loci of gaps in CGC2 after scaffolding. This is Supplementary Table 4.
- `dontFlip.ids`
	- CGC2 chromosome IDs that should not be reoriented.
- `flip.ids`
	- CGC2 chromosome IDs that need to be reoriented to match the orientation of other *C. briggsae* reference genomes.
- `HiC_contig_coverage.tsv`
	- Hi-C contact coverage of CGC2 contigs.
- `HiC_PB420_alignment.markedDUP.bam` This is a large file not present in this repo.
	- Alignment file (BAM) of Hi-C reads to the contiguous CGC2 genome with duplicate Hi-C reads marked. This file is created from `../../scripts/scaffolding/picard.sh`
- `HiC_PB420_alignment.sorted.bam` This is a large file not present in this repo.
	- Alignment file (BAM) of Hi-C reads to the contiguous CGC2 genome. This file is created from `../../scripts/scaffolding/bwa_hic_alignment.sh`
- `hifi_CGC2_noGaps_chroms_reoriented.bam` This is a large file not present in this repo.
	- HiFi read alignments to the final CGC2 genome. This file is created from `../../scripts/ont_hifi_CGC2_aln.sh`
- `hifi_CGC2_noGaps_chroms_reoriented.bam.bai`
	- Index file of associated BAM.
- `hifi_chromV_1kb_cov.regions.mosdepth.bed`
	- HiFi coverage per kb bin of CGC2 chromosome V. 
- `hifi_to_CGC2_scaffolds.bam` This is a large file not present in this repo.
	- HiFi read alignments to the scaffolded CGC2 genome with gaps and mis-oriented contigs present. This file is created from `../../scripts/ont_hifi_CGC2_aln.sh`
- `hifi_to_CGC2_scaffolds.bam.bai`
	- Index file of associated BAM.
- `ONT_CGC2_scaffolds.bam` This is a large file not present in this repo.
	- ONT read alignments to the scaffolded CGC2 genome with gaps mis-ordered contigs present. This file is created from `../../scripts/ont_hifi_CGC2_aln.sh`
- `ONT_CGC2_scaffolds.bam.bai`
	- Index file of associated BAM.
- `ONT_CGC2_noGap_chroms_reoriented.bam` This is a large file not present in this repo.
	- ONT read alignments to the final CGC2 genome. This file is created from `../../scripts/ont_hifi_CGC2_aln.sh`
- `ONT_CGC2_noGap_chroms_reoriented.bam.bai`
	- Index file of associated BAM.
- `ONT_chromV_1kb_cov.regions.mosdepth.bed`
	- ONT coverage per kb bin of CGC2 chromosome V.
- `out.hic`
	- Hi-C contact file output from Juicer.
- `out_JBAT.assembly`
	- Contiguous CGC2 assembly in a Juicebox-compatible format.
- `out_JBAT.hic`
	- Hi-C contact file otuput from Juice in a Juicebox-compatible format.
- `scaffold_contig_IDs.tsv`
	- Contig IDs and the corresponding scaffold IDs of where they were placed.
