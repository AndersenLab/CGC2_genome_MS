#!/bin/bash

#SBATCH -J AF16vcfFilt                    # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 8                           # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/biallelic_mito_nucl.oe  # Output log file
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/bialleic_mito.nucl.rr 

vcf="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/mutants.hard-filter.biallelic.NoMt.vcf.gz"
vcf_unfilt="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/mutants.hard-filter.vcf.gz"

bcftools view -m2 -M2 -e 'CHROM=="MtDNA"' -O z -o $vcf $vcf_unfilt
bcftools view -m2 -M2 -i 'CHROM=="MtDNA"' -O z -o ${vcf_unfilt%.vcf.gz}.biallelic.onlyMt.vcf.gz $vcf_unfilt 
