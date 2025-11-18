#!/bin/bash

#SBATCH -J snpeff_flatfile
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/slurm_output/snpeff_flatfile.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/slurm_output/snpeff_flatfile.rr 

SnpEff_annotated_vcf="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SnpEff/mutants.hard-filter.biallelic.onlyMt.snpeff.vcf.gz"
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/SnpEff"


if [[ ! -f  "$output_dir/SnpEff_flatFile.tsv" ]]; then 

    echo -e "Chromosome\tPosition\tREF\tALT\tconsequence\timpact\tAA_change\tALT_samples\ttranscript\tbackground_variant\tpossible_EMS" > $output_dir/SnpEff_flatFile.tsv

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN[\t%SAMPLE=%GT]\n' $SnpEff_annotated_vcf | \
    awk -F'\t' '{
        ALT_samples = "";  # Initialize string to collect samples with alt allele
        for (i = 6; i <= NF; i++) {  # Loop through fields containing sample=genotype
            if ($i ~ /0\/1|1\/0|1\/1/) {  # Check if sample has alt allele
                sub(/=.*/, "", $i);  # Remove genotype - leaving only sample name
                ALT_samples = ALT_samples (ALT_samples ? " " : "") $i;  # Append sample name
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
            split($5, annotations, ",");  # Split multiple annotations into array
            for (j in annotations) {
                split(annotations[j], snpEff, "|");  # Split each annotation by pipe

                SnpEff_consequence = snpEff[2];
                SnpEff_impact = snpEff[3];
                SnpEff_AA_change = (snpEff[11] != "" ? snpEff[11] : "N/A");
                transcript = (snpEff[4] != "" ? snpEff[4] : "N/A");
                gsub(/transcript:/, "", transcript); # Remove "transcript:" from the transcript string
                if (SnpEff_consequence == "intergenic_region") {
                    transcript = "N/A";
                }

                print $1"\t"$2"\t"$3"\t"$4"\t"SnpEff_consequence"\t"SnpEff_impact"\t"SnpEff_AA_change"\t"ALT_samples"\t"transcript"\t"background"\t"ems; 
            }
        }
    }' >> $output_dir/SnpEff_flatFile.tsv
else 
    echo "SnpEff flat file already exists."
fi


## adding two additional columns if a variant is found in every strain in relation to QX1410 (AF16 genetic background) and if the variant is potentially induced by EMS

