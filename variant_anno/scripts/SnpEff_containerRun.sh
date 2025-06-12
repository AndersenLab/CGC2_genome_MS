#!/bin/bash

#SBATCH -J snpEff                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 8:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 8                            # Number of cores
#SBATCH --mail-user=loconn13@jh.edu     # Email for job notifications
#SBATCH --mail-type=END                 # Notify when job ends
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/SnpEff_final.oe 
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/SnpEff_final.rr 

module load singularity  # Load Singularity module

# Define paths for binding directories
snpeff_input="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/raw_data/SnpEff"
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SnpEff"
containerImage="/vast/eande106/projects/Lance/THESIS_WORK/variant_annotation/container_images/loconn13999-snpeff_annotation_2024_09_19.sif"
vcf_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data"

# Use singularity exec for non-interactive execution
singularity exec --bind $snpeff_input:/snpeff_input \
                 --bind $output_dir:/annotation_output \
                 --bind $vcf_dir:/vcf_dir \
                 $containerImage /usr/bin/bash /vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/SnpEff_analysis.sh
