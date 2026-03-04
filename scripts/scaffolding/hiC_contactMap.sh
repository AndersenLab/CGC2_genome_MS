#!/bin/bash

bin="../../processed_data/genomes/PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed.bin"
agp="../../processed_data/genomes/PB420.20260223.inbred.withONT.blobFiltered.yahs_scaffoldeed_scaffolds_final.agp"
contigs="../../processed_data/genomes/PB420.20251025.inbred.withONT.blobFiltered.fa.fai"
scaff_sizes="../../processed_data/genomes/scaffold_sizes.tsv"

out="../../processed_data/scaffolding"

# Converting HiC coverage to a format compatible with Juicer and Juicebox using the juicer scripts from YaHS install
if [[ ! -s $out/alignments_sorted.txt || ! -s $out/out_JBAT.log ]]; then 
	$yahs_juicer pre $bin $agp $contigs \
		| sort -k2,2d -k6,6d -T $temp --parallel=24 -S32G \
		| awk 'NF' > $out/alignments_sorted.txt.part
	
	mv $out/alignments_sorted.txt.part $out/alignments_sorted.txt
	
	$yahs_juicer pre -a -o $out/out_JBAT $bin $agp $contigs > $out/out_JBAT.log 2>&1
else
	echo "First step complete"
fi

# Generate Hi-C contact map and file for Juicebox
if [[ ! -s $out/out.hic || ! -s $out/out_JBAT.hic ]]; then
	java -Xmx96G -jar $juicer pre $out/alignments_sorted.txt $out/out.hic.part $scaff_sizes 
	
	mv $out/out.hic.part $out/out.hic

	java -Xmx96G -jar $juicer pre $out/out_JBAT.txt $out/out_JBAT.hic.part <(cat $out/out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')
	mv $out/out_JBAT.hic.part $out/out_JBAT.hic
else
	echo "Hi-C contact map already created"
fi
