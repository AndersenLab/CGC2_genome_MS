#!/bin/bash

#SBATCH -J csq_flatfile
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/slurm_output/csq_flatfile.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/scripts/slurm_output/csq_flatfile.rr 

CSQ_annotated_vcf="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/CSQ/mutants.hard-filter.biallelic.NoMt.csq.vcf.gz"
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_genome_MS/variant_anno/processed_data/CSQ"
gff="/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/gff/c_briggsae.QX1410_20250929.csq.gff"

if [[ ! -f  "$output_dir/CSQ_flatFile.tsv" ]]; then 

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%SAMPLE=%GT]\n' $CSQ_annotated_vcf | \
    awk -F'\t' '{s
        ALT_samples = "";  
        for (i = 6; i <= NF; i++) { 
            if ($i ~ /0\/1|1\/0|1\/1/) { 
                sub(/=.*/, "", $i);  
                ALT_samples = ALT_samples (ALT_samples ? " " : "") $i; 
            }
        }
        
        n_ALT = split(ALT_samples, arr, " ")

        # print "ALT_samples length(chars)=" length(ALT_samples) ", count=" n_ALT > "/dev/stderr"
        
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

                    print $1"\t"$2"\t"$3"\t"$4"\t"CSQ_consequence"\t"CSQ_AA_change"\t"DNA_change"\t"ALT_samples"\t"transcript"\t"background"\t"ems;
                }
            } else {
                print $1"\t"$2"\t"$3"\t"$4"\tN/A\tN/A\tN/A\t"ALT_samples"\tN/A\t"background"\t"ems;
            }
        }
    }' >> $output_dir/CSQ_flatFile.tsv
else
    echo "CSQ flat file already exists."
fi



# Extracting transcripts and WBGeneID's from GFF3
# ID=transcript:Y74C9A.3.1;Parent=gene:WBGene00022277;Name=Y74C9A.3.1;wormpep=CE28146;locus=homt-1;uniprot_id=Q9N4D9;biotype=protein_coding
# Ensure ncRNAs are extracted as some annotation tools annotate these (e.g., VEP)
awk -F'\t' '$3 != "gene" && $3 != "intron" && $3 != "exon" && $3 != "five_prime_UTR" && $3 != "three_prime_UTR" && $3 != "CDS" && $3 != "start_codon" && $3 != "stop_codon" {print $9}' $gff | \
awk -F'[;=]' '{
    transcript=""; wbgene="";gene=""
    for(i=1; i<=NF; i++) {
        if ($i == "ID" && $(i+1) ~ /^transcript:/) { transcript=$(i+1); gsub("transcript:", "", transcript); }
        if ($i == "ID" && $(i+1) ~ /^Pseudogene:/) { transcript=$(i+1); gsub("Pseudogene:", "", transcript); }
        if ($i == "Parent" && $(i+1) ~ /^gene:/) { wbgene=$(i+1); gsub("gene:", "", wbgene); }
        if ($i == "sequence_name") { gene=$(i+1); }
        #if ($i == "Name") { name=$(i+1); } 
    }
    if (gene == "") { gene = "N/A"; }

    if (transcript != "" && wbgene != "") print transcript "\t" wbgene "\t" gene;
}' > $output_dir/$(basename $gff).tsv
sort -t$'\t' -k1,1 $output_dir/$(basename $gff).tsv > $output_dir/$(basename $gff)_sorted.tsv
rm $output_dir/$(basename $gff).tsv

# Adding WBGeneIDs and gene names
join -t$'\t' -1 1 -2 9 -a 2 $output_dir/$(basename $gff)_sorted.tsv <(sort -t$'\t' -k9,9 $output_dir/CSQ_flatFile.tsv) | \
        awk -F'\t' -v OFS='\t' '{
                if (NF == 13) {
                        print $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $1, $2, $3
                } else {
                        print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $1, $1  
                 }
        }' > $output_dir/CSQ_joined.tsv

sort -k1,1 -k2,2n $output_dir/CSQ_joined.tsv > $output_dir/CSQ_flatFile_WBGeneID.tsv
rm $output_dir/CSQ_flatFile.tsv $output_dir/CSQ_joined.tsv 






#### ADDING GRANTHAM SCORES ####
# CSQ 
awk -F'\t' -v OFS='\t' '
BEGIN {
    grantham["A/R"]=112; grantham["A/N"]=111; grantham["A/D"]=126; grantham["A/C"]=195; grantham["A/Q"]=91;
    grantham["A/E"]=107; grantham["A/G"]=60;  grantham["A/H"]=86;  grantham["A/I"]=94;  grantham["A/L"]=96;
    grantham["A/K"]=106; grantham["A/M"]=84;  grantham["A/F"]=113; grantham["A/P"]=27;  grantham["A/S"]=99;
    grantham["A/T"]=58;  grantham["A/W"]=148; grantham["A/Y"]=112; grantham["A/V"]=64;

    grantham["R/N"]=86;  grantham["R/D"]=96;  grantham["R/C"]=180; grantham["R/Q"]=43;  grantham["R/E"]=54;
    grantham["R/G"]=125; grantham["R/H"]=29;  grantham["R/I"]=97;  grantham["R/L"]=102; grantham["R/K"]=26;
    grantham["R/M"]=91;  grantham["R/F"]=97;  grantham["R/P"]=103; grantham["R/S"]=110; grantham["R/T"]=71;
    grantham["R/W"]=101; grantham["R/Y"]=77;  grantham["R/V"]=96;

    grantham["N/D"]=23;  grantham["N/C"]=139; grantham["N/Q"]=46;  grantham["N/E"]=42;  grantham["N/G"]=80;
    grantham["N/H"]=68;  grantham["N/I"]=149; grantham["N/L"]=153; grantham["N/K"]=94;  grantham["N/M"]=142;
    grantham["N/F"]=158; grantham["N/P"]=91;  grantham["N/S"]=46;  grantham["N/T"]=65;  grantham["N/W"]=174;
    grantham["N/Y"]=143; grantham["N/V"]=133;

    grantham["D/C"]=154; grantham["D/Q"]=61;  grantham["D/E"]=45;  grantham["D/G"]=94;  grantham["D/H"]=81;
    grantham["D/I"]=168; grantham["D/L"]=172; grantham["D/K"]=101; grantham["D/M"]=160; grantham["D/F"]=177;
    grantham["D/P"]=108; grantham["D/S"]=65;  grantham["D/T"]=85;  grantham["D/W"]=181; grantham["D/Y"]=160;
    grantham["D/V"]=152;

    grantham["C/Q"]=154; grantham["C/E"]=158; grantham["C/G"]=159; grantham["C/H"]=174; grantham["C/I"]=198;
    grantham["C/L"]=198; grantham["C/K"]=202; grantham["C/M"]=196; grantham["C/F"]=205; grantham["C/P"]=169;
    grantham["C/S"]=112; grantham["C/T"]=149; grantham["C/W"]=215; grantham["C/Y"]=194; grantham["C/V"]=192;

    grantham["Q/E"]=29;  grantham["Q/G"]=87;  grantham["Q/H"]=24;  grantham["Q/I"]=109; grantham["Q/L"]=113;
    grantham["Q/K"]=53;  grantham["Q/M"]=101; grantham["Q/F"]=116; grantham["Q/P"]=76;  grantham["Q/S"]=68;
    grantham["Q/T"]=42;  grantham["Q/W"]=130; grantham["Q/Y"]=99;  grantham["Q/V"]=96;

    grantham["E/G"]=98;  grantham["E/H"]=40;  grantham["E/I"]=134; grantham["E/L"]=138; grantham["E/K"]=56;
    grantham["E/M"]=126; grantham["E/F"]=140; grantham["E/P"]=93;  grantham["E/S"]=80;  grantham["E/T"]=65;
    grantham["E/W"]=152; grantham["E/Y"]=122; grantham["E/V"]=121;

    grantham["G/H"]=98;  grantham["G/I"]=135; grantham["G/L"]=138; grantham["G/K"]=127; grantham["G/M"]=127;
    grantham["G/F"]=153; grantham["G/P"]=42;  grantham["G/S"]=56;  grantham["G/T"]=59;  grantham["G/W"]=184;
    grantham["G/Y"]=147; grantham["G/V"]=109;

    grantham["H/I"]=94;  grantham["H/L"]=99;  grantham["H/K"]=32;  grantham["H/M"]=87;  grantham["H/F"]=100;
    grantham["H/P"]=77;  grantham["H/S"]=89;  grantham["H/T"]=47;  grantham["H/W"]=115; grantham["H/Y"]=83;
    grantham["H/V"]=84;

    grantham["I/L"]=5;   grantham["I/K"]=102; grantham["I/M"]=10;  grantham["I/F"]=21;  grantham["I/P"]=95;
    grantham["I/S"]=142; grantham["I/T"]=89;  grantham["I/W"]=61;  grantham["I/Y"]=33;  grantham["I/V"]=29;

    grantham["L/K"]=107; grantham["L/M"]=15;  grantham["L/F"]=22;  grantham["L/P"]=98;  grantham["L/S"]=145;
    grantham["L/T"]=92;  grantham["L/W"]=61;  grantham["L/Y"]=36;  grantham["L/V"]=32;

    grantham["K/M"]=95;  grantham["K/F"]=102; grantham["K/P"]=103; grantham["K/S"]=121; grantham["K/T"]=78;
    grantham["K/W"]=110; grantham["K/Y"]=85;  grantham["K/V"]=97;

    grantham["M/F"]=28;  grantham["M/P"]=87;  grantham["M/S"]=135; grantham["M/T"]=81;  grantham["M/W"]=67;
    grantham["M/Y"]=36;  grantham["M/V"]=21;

    grantham["F/P"]=114; grantham["F/S"]=155; grantham["F/T"]=103; grantham["F/W"]=40;  grantham["F/Y"]=22;
    grantham["F/V"]=50;

    grantham["P/S"]=74;  grantham["P/T"]=38;  grantham["P/W"]=147; grantham["P/Y"]=110; grantham["P/V"]=68;

    grantham["S/T"]=58;  grantham["S/W"]=177; grantham["S/Y"]=144; grantham["S/V"]=124;

    grantham["T/W"]=128; grantham["T/Y"]=92;  grantham["T/V"]=69;

    grantham["W/Y"]=37;  grantham["W/V"]=88;

    grantham["Y/V"]=55;
}

{
    if ($5 ~ /missense/) {
        
        # Extract amino acid substitution from column 6: e.g. 175D>175N -> D/N
        match($6, /[0-9]+([A-Z])>[0-9]*([A-Z])/, aa)
        ref = aa[1]
        alt = aa[2]

        key = ref "/" alt
        rev_key = alt "/" ref  ### Grantham matrix is symmetrical

        if (key in grantham)
            print $0, grantham[key]
        else if (rev_key in grantham)
            print $0, grantham[rev_key]
        else
            print $0, "N/A"
    } else {
        print $0, "N/A"
    }
}' $output_dir/CSQ_flatFile_WBGeneID.tsv > $output_dir/CSQ_WBGeneID_GRANTHAM_test.tsv




#### ADDING BLOSUM62 SCORES ####
# CSQ
awk -F'\t' -v OFS='\t' '
BEGIN {
    # BLOSUM62 substitution matrix 
    blosum["A/A"]=4;  blosum["A/R"]=-1; blosum["A/N"]=-2; blosum["A/D"]=-2; blosum["A/C"]=0;
    blosum["A/Q"]=-1; blosum["A/E"]=-1; blosum["A/G"]=0;  blosum["A/H"]=-2; blosum["A/I"]=-1;
    blosum["A/L"]=-1; blosum["A/K"]=-1; blosum["A/M"]=-1; blosum["A/F"]=-2; blosum["A/P"]=-1;
    blosum["A/S"]=1;  blosum["A/T"]=0;  blosum["A/W"]=-3; blosum["A/Y"]=-2; blosum["A/V"]=0;

    blosum["R/R"]=5;  blosum["R/N"]=0;  blosum["R/D"]=-2; blosum["R/C"]=-3; blosum["R/Q"]=1;
    blosum["R/E"]=0;  blosum["R/G"]=-2; blosum["R/H"]=0;  blosum["R/I"]=-3; blosum["R/L"]=-2;
    blosum["R/K"]=2;  blosum["R/M"]=-1; blosum["R/F"]=-3; blosum["R/P"]=-2; blosum["R/S"]=-1;
    blosum["R/T"]=-1; blosum["R/W"]=-3; blosum["R/Y"]=-2; blosum["R/V"]=-3;

    blosum["N/N"]=6;  blosum["N/D"]=1;  blosum["N/C"]=-3; blosum["N/Q"]=0;  blosum["N/E"]=0;
    blosum["N/G"]=0;  blosum["N/H"]=1;  blosum["N/I"]=-3; blosum["N/L"]=-3; blosum["N/K"]=0;
    blosum["N/M"]=-2; blosum["N/F"]=-3; blosum["N/P"]=-2; blosum["N/S"]=1;  blosum["N/T"]=0;
    blosum["N/W"]=-4; blosum["N/Y"]=-2; blosum["N/V"]=-3;

    blosum["D/D"]=6;  blosum["D/C"]=-3; blosum["D/Q"]=0;  blosum["D/E"]=2;  blosum["D/G"]=-1;
    blosum["D/H"]=-1; blosum["D/I"]=-3; blosum["D/L"]=-4; blosum["D/K"]=-1; blosum["D/M"]=-3;
    blosum["D/F"]=-3; blosum["D/P"]=-1; blosum["D/S"]=0;  blosum["D/T"]=-1; blosum["D/W"]=-4;
    blosum["D/Y"]=-3; blosum["D/V"]=-3;

    blosum["C/C"]=9;  blosum["C/Q"]=-3; blosum["C/E"]=-4; blosum["C/G"]=-3; blosum["C/H"]=-3;
    blosum["C/I"]=-1; blosum["C/L"]=-1; blosum["C/K"]=-3; blosum["C/M"]=-1; blosum["C/F"]=-2;
    blosum["C/P"]=-3; blosum["C/S"]=-1; blosum["C/T"]=-1; blosum["C/W"]=-2; blosum["C/Y"]=-2;
    blosum["C/V"]=-1;

    blosum["Q/Q"]=5;  blosum["Q/E"]=2;  blosum["Q/G"]=-2; blosum["Q/H"]=0;  blosum["Q/I"]=-3;
    blosum["Q/L"]=-2; blosum["Q/K"]=1;  blosum["Q/M"]=0;  blosum["Q/F"]=-3; blosum["Q/P"]=-1;
    blosum["Q/S"]=0;  blosum["Q/T"]=-1; blosum["Q/W"]=-2; blosum["Q/Y"]=-1; blosum["Q/V"]=-2;

    blosum["E/E"]=5;  blosum["E/G"]=-2; blosum["E/H"]=0;  blosum["E/I"]=-3; blosum["E/L"]=-3;
    blosum["E/K"]=1;  blosum["E/M"]=-2; blosum["E/F"]=-3; blosum["E/P"]=-1; blosum["E/S"]=0;
    blosum["E/T"]=-1; blosum["E/W"]=-3; blosum["E/Y"]=-2; blosum["E/V"]=-2;

    blosum["G/G"]=6;  blosum["G/H"]=-2; blosum["G/I"]=-4; blosum["G/L"]=-4; blosum["G/K"]=-2;
    blosum["G/M"]=-3; blosum["G/F"]=-3; blosum["G/P"]=-2; blosum["G/S"]=0;  blosum["G/T"]=-2;
    blosum["G/W"]=-2; blosum["G/Y"]=-3; blosum["G/V"]=-3;

    blosum["H/H"]=8;  blosum["H/I"]=-3; blosum["H/L"]=-3; blosum["H/K"]=-1; blosum["H/M"]=-2;
    blosum["H/F"]=-1; blosum["H/P"]=-2; blosum["H/S"]=-1; blosum["H/T"]=-2; blosum["H/W"]=-2;
    blosum["H/Y"]=2;  blosum["H/V"]=-3;

    blosum["I/I"]=4;  blosum["I/L"]=2;  blosum["I/K"]=-3; blosum["I/M"]=1;  blosum["I/F"]=0;
    blosum["I/P"]=-3; blosum["I/S"]=-2; blosum["I/T"]=-1; blosum["I/W"]=-3; blosum["I/Y"]=-1;
    blosum["I/V"]=3;

    blosum["L/L"]=4;  blosum["L/K"]=-2; blosum["L/M"]=2;  blosum["L/F"]=0;  blosum["L/P"]=-3;
    blosum["L/S"]=-2; blosum["L/T"]=-1; blosum["L/W"]=-2; blosum["L/Y"]=-1; blosum["L/V"]=1;

    blosum["K/K"]=5;  blosum["K/M"]=-1; blosum["K/F"]=-3; blosum["K/P"]=-1; blosum["K/S"]=0;
    blosum["K/T"]=-1; blosum["K/W"]=-3; blosum["K/Y"]=-2; blosum["K/V"]=-2;

    blosum["M/M"]=5;  blosum["M/F"]=0;  blosum["M/P"]=-2; blosum["M/S"]=-1; blosum["M/T"]=-1;
    blosum["M/W"]=-1; blosum["M/Y"]=-1; blosum["M/V"]=1;

    blosum["F/F"]=6;  blosum["F/P"]=-4; blosum["F/S"]=-2; blosum["F/T"]=-2; blosum["F/W"]=1;
    blosum["F/Y"]=3;  blosum["F/V"]=-1;

    blosum["P/P"]=7;  blosum["P/S"]=-1; blosum["P/T"]=-1; blosum["P/W"]=-4; blosum["P/Y"]=-3;
    blosum["P/V"]=-2;

    blosum["S/S"]=4;  blosum["S/T"]=1;  blosum["S/W"]=-3; blosum["S/Y"]=-2; blosum["S/V"]=-2;

    blosum["T/T"]=5;  blosum["T/W"]=-2; blosum["T/Y"]=-2; blosum["T/V"]=0;

    blosum["W/W"]=11; blosum["W/Y"]=2;  blosum["W/V"]=-3;

    blosum["Y/Y"]=7;  blosum["Y/V"]=-1;

    blosum["V/V"]=4;
}
{
    if ($5 ~ /missense/) {
        match($6, /[0-9]+([A-Z])>[0-9]*([A-Z])/, aa)
        ref = aa[1]
        alt = aa[2]

        key = ref "/" alt
        rev_key = alt "/" ref # BLOSUM62 matrix is symmetrical

        if (key in blosum)
            print $0, blosum[key]
        else if (rev_key in blosum)
            print $0, blosum[rev_key]
        else
            print $0, "N/A"
    } else {
        print $0, "N/A"
    }
}' $output_dir/CSQ_WBGeneID_GRANTHAM_test.tsv > $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_test.tsv






# ### ADDING PERCENT PROTEIN #### 
# ----------- Load CDS info into memory ------------
declare -A CDS_map

while IFS=$'\t' read -r chr type start end strand attr; do
    transcript=$(echo "$attr" | grep -o 'Parent=transcript:[^;]*' | cut -d: -f2)
    # echo "Processing transcript: $transcript - Adding to CDS_map associative array"
    CDS_map["$transcript"]+="$chr,$start,$end,$strand;"
done < <(awk -F'\t' -v OFS='\t' '$3 == "CDS" {print $1, $3, $4, $5, $7, $9}' $gff)

echo "Made associative array for ${#CDS_map[@]} CDS entries"


### CSQ 
#Create header only if the output doesn't already exist
if [[ ! -f $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_PP_final.tsv ]]; then
    # echo "chrom,pos,ref,alt,consequence,AA,DNAchange,strain,divergent,transcript_name,wbgene,gene_name,grantham_score,blosum_score,percent_protein" > $output_dir/CSQ_${1}_WBGeneID_GRANTHAM_BLOSUM_PP_final.csv
        echo -e "Chromosome\tPosition\tREF\tALT\tconsequence\tAA_change\tDNA_change\tALT_samples\ttranscript\tbackground_variant\tpossible_EMS\twbgene\tgene_name\tgrantham_score\tblosum_score\tpercent_protein" > $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_PP_final.tsv
else 
    echo "Final CSQ file for PP calculation is already created."
fi

# Create an associative array of processed keys for quick lookup
declare -A processed

# Populate the processed set from existing output
tail -n +2 $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_PP_final.tsv | awk -F',' '{print $2"_"$11}' | while read key; do 
    processed["$key"]=1
done

while IFS=$'\t' read -r chrom pos ref alt consequence AA DNAchange strain background ems transcript_name wbgene locus grantham_score blosum_score; do  

    key="${pos}_${transcript_name}"
    if [[ -n "${processed[$key]}" ]]; then
        echo "PP already calculated for ${processed[$key], continuing...}"
        continue
    fi

    pp="N/A"

    if [[ $consequence =~ missense|frameshift|stop|deletion|synonymous|insertion ]]; then
        cds_entries="${CDS_map[$transcript_name]}"
        # echo "Consequence is: $consequence. Calculating PP" # this is working correctly
        
        cds_entries="${CDS_map[$transcript_name]}"
        if [[ -n "$cds_entries" ]]; then
            # echo "Consequence is: $consequence. Calculating PP for $cds_entries"

            IFS=';' read -ra cds_array <<< "$cds_entries"
            found=0
            cum_len=0
            CDS_total_len=0

            # Total CDS length first
            for entry in "${cds_array[@]}"; do
                IFS=',' read -r _ start end _ <<< "$entry"
                CDS_total_len=$(( CDS_total_len + (end - start + 1) ))
            done

            # Extract strand info
            IFS=',' read -r _ _ _ strand <<< "${cds_array[0]}"

            if [[ $strand == "-" ]]; then
                for ((i=${#cds_array[@]}-1; i>=0; i--)); do
                    IFS=',' read -r chr start end _ <<< "${cds_array[i]}"
                    cds_segment=$(( end - start + 1 ))

                    # echo "Comparing pos=$pos to start=$start and end=$end"

                    if (( pos >= start && pos <= end )); then
                        nucl_number=$(( cum_len + (end - pos + 1) ))
                        found=1
                        break
                    fi
                    cum_len=$(( cum_len + cds_segment ))
                done
            else
                for entry in "${cds_array[@]}"; do
                    IFS=',' read -r chr start end _ <<< "$entry"
                    cds_segment=$(( end - start + 1 ))

                    # echo "Comparing pos=$pos to start=$start and end=$end"

                    if (( pos >= start && pos <= end )); then
                        nucl_number=$(( cum_len + (pos - start + 1) ))
                        found=1
                        break
                    fi
                    cum_len=$(( cum_len + cds_segment ))
                done
            fi

            if (( found && CDS_total_len > 0 )); then
                pp=$(awk -v r=$nucl_number -v t=$CDS_total_len 'BEGIN { printf "%.2f", (r / t * 100) }')
            fi
        fi
    fi

    echo -e "$chrom\t$pos\t$ref\t$alt\t$consequence\t$AA\t$DNAchange\t$strain\t$background\t$ems\t$transcript_name\t$wbgene\t$locus\t$grantham_score\t$blosum_score\t$pp" >> $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_PP_final.tsv

done < $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_test.tsv

rm $output_dir/CSQ_flatFile_WBGeneID.tsv $output_dir/c_briggsae.QX1410_20250929.csq.gff_sorted.tsv $output_dir/CSQ_WBGeneID_GRANTHAM_BLOSUM_test.tsv $output_dir/CSQ_WBGeneID_GRANTHAM_test.tsv 






