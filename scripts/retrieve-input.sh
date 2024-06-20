#!/usr/bin/env bash

set -euo pipefail

input_dir="$1"

output_fastq="$2"
output_fasta="$3"

# TODO For older samples the data might be in a differently-named file, could
# branch to handle that here.
tar -x -f "$input_dir/basecalled_output.tar" \
    --to-stdout \
    --wildcards 'basecalled_output/pass/*.fastq.gz' \
| gunzip - \
> "$output_fastq"

# Example for the test dataset
# cat "$input_dir"/*.fastq > "$output_fastq"

awk '
    NR % 4 == 1 { printf(">%s\n",substr($0,2)) }
    NR % 4 == 2 { print }
' "$output_fastq" > "$output_fasta"
