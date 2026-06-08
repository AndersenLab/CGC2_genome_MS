#!/usr/bin/env bash

query="$1"
DB="$2"

blastn -task blastn-short \
  -query $query \
  -db $DB \
  -dust no \
  -soft_masking false \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand bitscore evalue"
