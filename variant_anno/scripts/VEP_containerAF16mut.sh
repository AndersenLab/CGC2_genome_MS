#!/bin/bash

#SBATCH -J VEP_troubleshoot       # Job name
#SBATCH -A eande106          # Allocation name
#SBATCH -p parallel          # Partition/Queue name
#SBATCH -t 12:00:00          # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                 # Number of nodes
#SBATCH -n 12                 # Number of cores
#SBATCH --mail-user=loconn13@jh.edu  # Email for notifications
#SBATCH --mail-type=END      # Notify when job ends
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/VEP.oe # Output log
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/VEP.rr  # Error log

# Container image and environment setup
container_image="/vast/eande106/projects/Lance/THESIS_WORK/variant_annotation/container_images/loconn13999-ensemble_vep_2024_05_24.sif"
raw_data="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/raw_data"
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/VEP"
vcf_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data"

gff="/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/gff/c_briggsae.QX1410_20250929.csq.gff"
ref_genome="c_briggsae.QX1410_nanopore.Feb2020.genome.fa"
vcf="mutants.hard-filter.biallelic.NoMt.vcf.gz"

module load singularity

mkdir -p $raw_data/VEP
cd $raw_data/VEP

# Sort GFF file for format that VEP needs
grep -v "#" $gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > "$(basename $gff .gff3)_VEPsorted.gff3.gz"
tabix -p gff "$(basename $gff .gff3)_VEPsorted.gff3.gz"

mkdir -p $output_dir

# Run with only local files, no cache and no downloads
singularity exec --bind "$raw_data/VEP:/raw_data" \
    --bind "$output_dir:/annotation_output" \
    --bind /vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020:/vast_data \
    --bind $vcf_dir:/vcf \
    $container_image vep \
    --vcf \
    --force_overwrite \
    --input_file /vcf/$vcf \
    --fasta /vast_data/$ref_genome \
    --gff /raw_data/$(basename $gff .gff3)_VEPsorted.gff3.gz \
    --output_file /annotation_output/$(basename $vcf .vcf.gz).VEP.vcf.gz \
    --compress_output bgzip 
