#!/usr/bin/env bash

set -euo pipefail

out_blast="$1"
shift

out_cigar="$1"
shift

bam="$1"
shift

# Remaining arguments are the blast results for the various chromosomes.

# Collect those into a single file here.
cat "$@" > "$out_blast"

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
