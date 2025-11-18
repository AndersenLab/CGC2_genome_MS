#!/bin/bash

#SBATCH -J csq_flatfile
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/variant_annotation/processed_data/flat_file_creation/c_briggsae/SLURM_output/0729_VEP.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/variant_annotation/processed_data/flat_file_creation/c_briggsae/SLURM_output/0729_VEP.rr 

CSQ_annotated_vcf="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/CSQ/mutants.hard-filter.biallelic.NoMt.csq.vcf.gz"
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/CSQ"

if [[ ! -f  "$output_dir/CSQ_flatFile.tsv" ]]; then 
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%SAMPLE=%GT]\n' $CSQ_annotated_vcf | \
    awk -F'\t' '{
        ALT_samples = "";  
        for (i = 6; i <= NF; i++) { 
            if ($i ~ /0\/1|1\/0|1\/1/) { 
                sub(/=.*/, "", $i);  
                ALT_samples = ALT_samples (ALT_samples ? " " : "") $i; 
            }
        }

        if (ALT_samples != "") {  
            if ($5 != ".") {
                split($5, csq_entries, ",");  # Split multiple annotations by ","

                for (j in csq_entries) {
                    split(csq_entries[j], csq, "|");  # Extract annotation components
                    
                    if (csq[1] ~ /^@/) {
                        CSQ_consequence = csq[1];  
                        CSQ_AA_change = csq[1]; 
                        DNA_change = csq[1];
                        transcript = csq[1];  
                    } else {
                        CSQ_consequence = (csq[1] != "" ? csq[1] : "N/A");
                        CSQ_AA_change = (csq[6] != "" ? csq[6] : "N/A");
                        DNA_change = (csq[7] != "" ? csq[7] : "N/A");
                        transcript = (csq[3] != "" ? csq[3] : "N/A");
                    }

                    print $1"\t"$2"\t"$3"\t"$4"\t"CSQ_consequence"\t"CSQ_AA_change"\t"DNA_change"\t"ALT_samples"\t"transcript;
                }
            } else {
                print $1"\t"$2"\t"$3"\t"$4"\tN/A\tN/A\tN/A\t"ALT_samples"\tN/A";
            }
        }
    }' >> $output_dir/CSQ_flatFile.tsv
else
    echo "CSQ flat file already exists."
fi


## adding two additional columns if a variant is found in every strain in relation to QX1410 (AF16 genetic background) and if the variant is potentially induced by EMS


