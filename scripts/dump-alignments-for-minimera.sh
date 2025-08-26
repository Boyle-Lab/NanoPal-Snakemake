#!/usr/bin/env bash

set -euo pipefail

bam="$1"

# Flags:
# 0x004 Exclude unmapped reads
# 0x100 Exclude secondary alignments
# 0x400 Exclude PCR/optical duplicates
# 0x800 Exclude supplementary alignments

samtools view "$bam" \
    -F 0x004 \
    -F 0x100 \
    -F 0x400 \
    -F 0x800 \
| gawk '
    BEGIN {
        OFS=","
        print "id", "contig", "pos", "rev", "flag", "cigar"
    }

    {
        id     = $1
        flag   = $2
        contig = $3
        pos    = $4
        cigar  = $6
        rev = rshift(and(flag, 0x10), 4) # 0x10 is the rev bit
        print id, contig, pos, rev, flag, cigar
    }
'
