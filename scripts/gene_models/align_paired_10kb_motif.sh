#submitted as sbatch job to HPC
#activate conda environment
source activate star


#$wkdir is working directory
#$prefix is sample ID / name
#$GENOME is the reference genome
#Aligned trimmed FASTQ to reference genome
mkdir -p $wkdir/processed_data/alignments_motif/$prefix/

cd $wkdir/processed_data/alignments_motif/$prefix/
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 600000000000 \
--genomeDir . \
--genomeFastaFiles $wkdir/raw_data/genome/$GENOME \
--genomeSAindexNbases 12 \
--alignIntronMax 10000 \
--outSAMstrandField intronMotif
STAR \
--runThreadN 24 \
--genomeDir . \
--outSAMtype BAM Unsorted SortedByCoordinate \
--twopassMode Basic \
--readFilesCommand zcat \
--alignIntronMax 10000 \
--outSAMstrandField intronMotif \
--readFilesIn $wkdir/processed_data/trim/$prefix/${prefix}_R1_001.fastq.trimmed.gz  $wkdir/processed_data/trim/$prefix/${prefix}_R2_001.fastq.trimmed.gz

