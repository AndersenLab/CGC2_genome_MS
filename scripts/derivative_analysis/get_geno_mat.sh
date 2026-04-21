bcftools query -f '%CHROM:%POS[\t%GT]\n' all_derivative.hard.vcf.gz > hard.gt_raw.tsv
bcftools query -l all_derivative.hard.vcf.gz > samples.txt
printf "Site\t" > header.txt
paste -sd '\t' samples.txt >> header.txt
printf "\n" >> header.txt
awk 'BEGIN{OFS="\t"}
{
  printf $1;
  for(i=2;i<=NF;i++){
    gt=$i;
    if(gt=="0/0" || gt=="0|0") val=0;
    else if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") val=1;
    else if(gt=="1/1" || gt=="1|1") val=2;
    else val="NA";
    printf OFS val;
  }
  printf "\n";
}' hard.gt_raw.tsv > hard.gt_numeric.tsv
cat header.txt hard.gt_numeric.tsv > genotype_matrix.tsv
