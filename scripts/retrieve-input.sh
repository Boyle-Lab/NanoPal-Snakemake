#!/usr/bin/env bash

set -euo pipefail

input_dirs="$1"

output_fastq="$2"
output_fasta="$3"
output_basecall_info="$4"

rm -f "$output_fastq"
rm -f "$output_basecall_info"

function retrieve_input {
    # TODO For older acquisitions the data might be in a differently-named file,
    # could branch to handle that here.
    tar="$1/basecalled_output.tar"

    # Pull out all the reads into a FASTQ.
    tar -x -f "$tar" \
        --to-stdout \
        --wildcards 'basecalled_output/pass/*.fastq.gz' \
    | gunzip - \
    >> "$output_fastq"

    # Pull out the basecalling metadata into a separate text file.
    log=$(tar --list -f "$tar" \
              --wildcards 'basecalled_output/ont_basecall_client_log-*' \
          | sort | head -1)
    tar -x -f "$tar" \
        --to-stdout \
        "$log" \
    | head -n 20 \
    >> "$output_basecall_info"

    # Example for the test dataset
    # cat "$input_dir"/*.fastq > "$output_fastq"
}

export -f retrieve_input

echo "$input_dirs" \
    | xargs -I DIR bash -c 'retrieve_input "$@"' _ DIR

# Convert the FASTQ to a FASTA.
awk '
    NR % 4 == 1 { printf(">%s\n",substr($0,2)) }
    NR % 4 == 2 { print }
' "$output_fastq" > "$output_fasta"
