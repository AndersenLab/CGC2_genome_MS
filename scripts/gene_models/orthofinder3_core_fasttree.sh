#submitted as sbatch job to HPC
#activate environment
source activate of3_env

#$tag is output results folder name
#$prot is the folder of protein sequences
#run orthofinder
orthofinder -M msa -S diamond_ultra_sens -A famsa -n $tag -t 48 -f $prot
