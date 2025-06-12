#!/bin/bash

#SBATCH -J CSQ_analysis                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 8:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 12                            # Number of cores per task
#SBATCH --mail-user=loconn13@jh.edu     # Email for job notifications
#SBATCH --mail-type=END                 # Notify when job ends
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/CSQ.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SLURM_output/CSQ.rr 

module load singularity  # Load Singularity module

# Define paths for binding directories
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/CSQ"
containerImage="/vast/eande106/projects/Lance/THESIS_WORK/variant_annotation/container_images/loconn13999-csq_annotation_2024_09_18.sif"
data_dir="/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020"

vcf="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/mutants.hard-filter.biallelic.NoMt.vcf.gz"
gff="$data_dir/csq/QX1410.update.April2025.noWBGeneID.csq.gff3"
ref_genome="$data_dir/c_briggsae.QX1410_nanopore.Feb2020.genome.fa"

mkdir -p $output_dir

# Use singularity exec for non-interactive execution
singularity exec --bind $output_dir:/annotation_output \
                 --bind $data_dir:/data_dir \
                 $containerImage bcftools csq \
                 -O z --fasta-ref $ref_genome \
                 --gff-annot $gff \
                 --ncsq 1000 \
		 --phase a $vcf > "$output_dir/$(basename ${vcf} .vcf.gz).csq.vcf.gz"

