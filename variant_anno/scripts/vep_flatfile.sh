#!/bin/bash

#SBATCH -J vep_flatfile
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/slurm_output/vep_flatfile.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/slurm_output/vep_flatfile.rr 

VEP_annotated_vcf="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/VEP/mutants.hard-filter.biallelic.NoMt.VEP.vcf.gz"
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/VEP"



if [[ ! -f  "$output_dir/VEP_flatFile.tsv" ]]; then 

    echo -e "Chromosome\tPosition\tREF\tALT\tconsequence\timpact\tAA_change\tALT_samples\ttranscript\tbackground_variant\tpossible_EMS" > $output_dir/VEP_flatFile.tsv
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ[\t%SAMPLE=%GT]\n' "$VEP_annotated_vcf" | \
    awk -F'\t' '{
        ALT_samples = "";  
        for (i = 6; i <= NF; i++) {  
            if ($i ~ /0\/1|1\/0|1\/1/) { 
                sub(/=.*/, "", $i);  
                ALT_samples = ALT_samples (ALT_samples ? " " : "") $i; 
            }
        }

        n_ALT = split(ALT_samples, arr, " ")
        
        ## adding two additional columns if a variant is found in every strain in relation to QX1410 (AF16 genetic background) and if the variant is potentially induced by EMS
        if (n_ALT == 20) {
            background = "yes"
        } else {
            background = "no"
        }

        if ($3 == "G" && $4 == "A" || $3 == "C" && $4 == "T") {
            ems = "yes"
        } else {
            ems = "no"
        }

        if (ALT_samples != "") {
            if ($5 != ".") {
                split($5, annotations, ",");
                for (j in annotations) {
                    split(annotations[j], vep, "|");  #split VEP annotation by pipe
                    VEP_consequence = vep[2];  
                    VEP_impact = vep[3];  
                    VEP_AA_change = (vep[16] != "" ? vep[16] : "N/A");
                    transcript = (vep[7] != "" ? vep[7] : "N/A");
                    if (VEP_consequence == "intergenic_variant") {
                        transcript = "N/A";
                    }
                    if (transcript ~ /&/) {
                        sub(/^.*&/, "", transcript);
                    }
                    # blosum62score = (vep[25] != "" ? vep[25] : "N/A");
                    print $1"\t"$2"\t"$3"\t"$4"\t"VEP_consequence"\t"VEP_impact"\t"VEP_AA_change"\t"ALT_samples"\t"transcript"\t"background"\t"ems; 
                }
            } else {
                print $1"\t"$2"\t"$3"\t"$4"\tN/A\tN/A\tN/A\t"ALT_samples"\tN/A"background"\t"ems;
            }
        }
    }' >> $output_dir/VEP_flatFile.tsv
else 
    echo "VEP flat file already exists."
fi


## adding two additional columns if a variant is found in every strain in relation to QX1410 (AF16 genetic background) and if the variant is potentially induced by EMS
