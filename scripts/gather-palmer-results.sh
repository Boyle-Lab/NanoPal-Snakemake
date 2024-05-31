#!/usr/bin/env bash

set -euo pipefail

palmer_dir="$1"
bam="$2"

out_blast="$3"
out_cigar="$4"

find "$palmer_dir" -name 'blastn_refine.txt' | xargs cat > "$out_blast"

samtools view "$bam" \
    --min-MQ 10 \
    --exclude-flags=0x100 \
    --exclude-flags=0x200 \
    --exclude-flags=0x400 \
    --exclude-flags=0x700 \
| awk '{{print $1, $3, $4, $6}}' > "$out_cigar"
#              read id, chromosome, position, cigar
