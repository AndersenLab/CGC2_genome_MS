#!/bin/bash

#SBATCH -J qc                    # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 2:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 6                           # Number of cores

output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data"

VEP_lines=$(cut -f1,2 $output_dir/VEP/VEP_flatFile.tsv | sort | uniq | wc -l)
CSQ_lines=$(cut -f1,2 $output_dir/CSQ/CQQ_flatFile.tsv | sort | uniq | wc -l)

if [[ $VEP_lines == $CSQ_lines ]]; then
    echo "Analysis complete: VEP and CSQ have the same number of variant calls"
else 
    echo "There is a mismatch in number of unique variant calls among the nuclear variant annotation tools: VEP=$VEP_lines; CSQ=$CSQ_lines."
