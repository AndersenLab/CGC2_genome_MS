## Scripts 
This directory contains the scripts that were used to perform analyses and generate figures and tables for the manuscript.

## Directory Structure

scripts/
- AF16cb5_QX1410_dotplot/
- AF16_CGC2_dotplot/
- AF16_QX1410_dotplot/
- CGC2_QX1410_dotplot/
- derivative_analysis/
- gene_models/
- marker_liftover/ 
- QX1410_QX1410_dotplot
- scaffolding/ 
- telomeres/ 
- VX34_QX1410_dotplot

## AF16cb5_QX1410_dotplot/
- `AF16cb5_QX1410_dotplot.R`
	- This script visualizes genome-genome alignments of AF16 cb5 to QX1410 and creates Supplementary Figure 2. 
- `AF16cb5_QX1410_nucmer.sh`  
	- Execution of nucmer to align the genome of AF16 cb5 to QX1410.

## AF16_CGC2_dotplot/
- `AF16_CGC2_dotplot.R`
	- This script visualizes genome-genome alignments of AF16 cb4 to CGC2 and creates Supplementary Figure 13.
- `AF16_CGC2_nucmer.sh`
	- Execution of nucmer to align the genome of AF16 cb4 to CGC2.

## AF16_QX1410_dotplot/
- `AF16_QX1410_dotplot.R`
	- This script visualizes genome-genome alignments of AF16 cb4 to QX1410 and creates Fig 1.
- `AF16_QX1410_nucmer.sh`
	- Execution of nucmer to align the genome of AF16 cb4 to QX1410. The output from this script is Supplementary Table 1.

## CGC2_QX1410_dotplot/
- `CGC2_QX1410_dotplot.R`
	- This script visualizes genome-genome alignments of CGC2 to QX1410 and creates Fig 2.
- `CGC2_QX1410_nucmer.sh`
	- Execution of nucmer to align the genome of CGC2 to QX1410. The output of this script is Supplementary Table 5.
- `CGC2_withGaps_QX1410_dotplot.R`
	- This script visualizes genome-genome alignments of CGC2 to the left-end of QX1410 chromosome V and creates Supplementary Fig 7.

## derivative_analysis/
- `get_geno_mat.sh`
	- Script that extracts genotype matrix from VCF
- `visualize_similarity.R`
	- Script that visualizes SNV differences between AF16 derivatives. Makes Supplementary Figure 1.

## gene_models/
- `align_paired_10kb_motif.sh`
	- Script to align RNAseq data with STAR.
- `fastp_trim.sh`
	- Script to trim RNAseq data.
- `orthofinder3_core_fasttree.sh`
	- Script to run orthofinder from protein sequences of CGC2, N2, QX1410, and AF16 cb4.
- `visualize_SCortho_3gen.R`
	- Script to visualize mapped and unmapped single copy orthologs between C. briggsae gene annotations.

## marker_liftover/
- `marker2fasta.R`
	- Script to convert marker TSV to primer sequence FASTA. 
- `process_blast.R`
	- Script to process BLAST results into mappable marker positions between AF16 cb4 and CGC2
- `blastn_short.sh`
	- Script to BLAST short primer sequences against a genome database.

## QX1410_QX1410_dotplot/
- `QX1410_QX1410_nucmer.sh`
	- Execution of nucmer to align the genome of QX1410 to itself. 
- `QX_QX_dotplot_INV.R`
	- Visualization of the QX1410 aligned to itself at the private inversion locus. The dotplot figure is displayed in Supplementary Figure 10a.

## scaffolding/
- `bwa_hic_alignment.sh`
	- Alignment of Hi-C seq data to contiguous CGC2 genome assembly. 
- `chromOrientation_flipping.sh`
	- The removal of the duplicate haplotig, scaffold 7, and flipping the orientation of CGC2 chromosomes to match the orientation of other *C. briggsae* reference genome chromosomes.
- `chrom_size_differences.R`
 	- Quantifying the chromosome size differences between CGC2 and AF16 cb4. This script creates Supplementary Figure 8.
- `coverage_rDNA_repeat.R`
	- This script visualizes the HiFi and ONT seq coverage over rDNA cistron repeat units found on the left-end of chromosome V in CGC2. This script creates Supplementary Figure 12.
- `duplicate_haplotig.R`
	- This script visualizes Hi-C seq coverage of each CGC2 contig and the alignment of the duplicate haplotig, scaffold 7, to another CGC2 scaffold. This script creates Supplementary Figure 4.
- `hiC_contactMap.sh`
	- This script creates the neccessary files for visualzing the CGC2 Hi-C contact map with Juicebox. The files created from this script are used for Hi-C contact map visualization displayed in Supplementary Figure 6.
- `mosdepth_coverage_rDNA.sh`
	- Quantification of HiFi and ONT seq coverage per kb bin of the CGC2 genome. 
- `ont_hifi_CGC2_aln.sh`
	- Alignment of HiFi and ONT seq data to CGC2 (pre and post gap closing) using minimap2.
- `picard.sh`
	- Marking duplicate Hi-C reads using picard for scaffolding with YaHS.
- `tgs_gapcloser.sh`
	- TGS-GapCloser to close the five gaps on the left-end of chromosome V in CGC2 using ONT reads and polishing with short-read data using Pilon.
- `yahs.sh`
	- Scaffolding of CGC2 using Hi-C reads with YaHS.

## telomeres/
- `telomere_density.R`
	- Visualization of telomeric repeat sequence in *C. briggsae* reference genomes. This script creates Supplementary Figure 11.
- `telomeres.sh`
	- Quantification of telomeric repeat sequence in *C. briggsae* reference genomes. 

## VX34_QX1410_dotplot/
- `VX34_CGC2_AF16_QX1410_chromV_inversion_vis.R`
	- Visualization of genome-genome alignments of *C. briggsae* reference genomes to QX1410 to visualize a ~625 kb inversion. This script creates Supplementary Figure 9.
- `VX34_QX1410_nucmer.sh`
	- Execution of nucmer to align the genome of VX34 to QX1410.
