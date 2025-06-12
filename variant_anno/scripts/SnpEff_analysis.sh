#!/bin/bash

vcf="/vcf_dir/mutants.hard-filter.biallelic.onlyMt.vcf.gz"
config="/snpeff_input/snpEff.config"
database_name="c_briggsae.QX1410_nanopore.Feb2020"

output_file="/annotation_output/$(basename ${vcf} .vcf.gz).snpeff.vcf.gz"

# Build the database
cd /snpeff_input #change to the data directory specified 

if [[ ! -f /snpeff_input/$database_name/snpEffectPredictor.bin ]]; then
    if ! java -jar /usr/bin/snpEff/snpEff/snpEff.jar build -noCheckCds -noCheckProtein -gtf22 -v $database_name; then
        echo "Error building SnpEff database for QX1410"
        exit 1
    fi
fi

# Run SnpEff annotation
if ! bcftools view -O v $vcf | \
   java -jar /usr/bin/snpEff/snpEff/snpEff.jar eff -csvStats /annotation_output/containerRun/snpeff.stats.csv \
              -nodownload \
              -dataDir /snpeff_input \
              -config $config \
              $database_name | \
   bcftools view -O z > $output_file; then
    echo "Error during VCF annotation"
    exit 1
fi
