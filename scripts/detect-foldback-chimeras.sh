#!/usr/bin/env bash

set -euo pipefail

threads="$1"
bam="$2"
result="$3"

function detect_chimera {
    set -euo pipefail

    id=$(echo "$1" | cut -d $'\x1f' -f 1)
    r=$(echo "$1"  | cut -d $'\x1f' -f 2)
    export r

    water \
        asis::"$r" \
        <(revseq -sequence asis::"$r" -outseq /dev/stdout) \
        -gapopen 5.0 -gapextend 0.5 \
        -outfile /dev/stdout \
        | grep Identity | tr -d '()%' | tr '/' ' ' | awk -v OFS=$'\t' -v id="$id" '{print id, $3, $4, $5}'
}

export -f detect_chimera

samtools fasta -F 0x100 -F 0x200 -F 0x400 -F 0x800 "$bam" \
    | paste -sd $'\x1f\n' \
    | xargs --max-procs "$threads" -n1 bash -c 'detect_chimera "$1"' - 2>/dev/null \
    > "$result"

