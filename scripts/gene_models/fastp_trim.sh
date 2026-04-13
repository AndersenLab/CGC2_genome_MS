#submitted as sbatch job to HPC
#activate conda environment
source activate fastp

#$wkdir is working directory
#$prefix is sample ID / name
#gets untrimmed FASTQ and returns trimmed FASTQ
mkdir -p $wkdir/processed_data/trim/$prefix/

fastp \
  -i $wkdir/raw_data/short_read/${prefix}_R1_001.fastq.gz \
  -I $wkdir/raw_data/short_read/${prefix}_R2_001.fastq.gz \
  -o $wkdir/processed_data/trim/$prefix/${prefix}_R1_001.fastq.trimmed.gz \
  -O $wkdir/processed_data/trim/$prefix/${prefix}_R2_001.fastq.trimmed.gz \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --trim_poly_x \
  --length_required 30 \
  --thread 24 \
  --html sample_fastp.html \
  --json sample_fastp.json
