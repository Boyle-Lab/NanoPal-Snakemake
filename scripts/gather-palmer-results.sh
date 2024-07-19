#!/usr/bin/env bash

set -euo pipefail

palmer_dir="$1"
bam="$2"

out_blast="$3"
out_cigar="$4"

# Grab all the Palmer blast results from the regional subsets in the working
# directories.
find "$palmer_dir" -name 'blastn_refine.txt' \
    | xargs cat > "$out_blast"

# Clean up Palmer directory once we've gotten what we need from it to avoid
# exhausting all the inodes on the drive when running many datasets.
find "$palmer_dir" -maxdepth 2 -name 'chr*_*_*' -type d \
    | xargs rm -r

# TODO Extract this to an earlier step, since it only depends on the alignment
# and it's wasteful to do it for every dataset+mei.
samtools view "$bam" \
    -q 10 \
    -F 0x100 \
    -F 0x200 \
    -F 0x400 \
    -F 0x800 \
| awk '{print $1, $3, $4, $6}' > "$out_cigar"
#        read id, chromosome, position, cigar
