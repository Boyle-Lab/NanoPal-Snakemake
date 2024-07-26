#!/usr/bin/env bash

set -euo pipefail
set -x

input_files="$1"

export output_fastq="$2"
export output_fasta="$3"
# export output_basecall_info="$4"

rm -f "$output_fastq"
# rm -f "$output_basecall_info"

function retrieve_from_tar {
    tar="$1"

    # Pull out all the reads into a FASTQ.
    tar -x -f "$tar" \
        --to-stdout \
        --wildcards 'basecalled_output/pass/*.fastq.gz' \
    | gunzip - \
    >> "$output_fastq"
}

function retrieve_from_bam {
    bam="$1"

    # Copy the tags from Dorado, see https://github.com/nanoporetech/dorado/blob/master/documentation/SAM.md
    samtools fastq "$bam" \
      -T MM,pi,sp,ns,ts,RG,qs,ts,ns,mx,ch,rn,st,du,fn,sm,sd,sv,mv,dx,pi,sp,pt,bh,MN \
    >> "$output_fastq"
}

function retrieve_input {
    set -euo pipefail
    set -x

    case "$1" in
      *.bam)
        retrieve_from_bam "$1"
        ;;
      *.tar | *.tar.gz | *.tgz)
        retrieve_from_tar "$1"
        ;;
      *)
        echo "Bad input file $1 -- expected BAM or tar." >&2
        exit 1
        ;;
    esac

    # TODO: Pull out the basecalling metadata into a separate text file.
    # log=$(tar --list -f "$tar" \
    #           --wildcards 'basecalled_output/ont_basecall_client_log-*' \
    #       | sort | head -n 1)

    # tar -x -f "$tar" --to-stdout "$log" | head -n 20 >> "$output_basecall_info"
}

export -f retrieve_input
export -f retrieve_from_bam
export -f retrieve_from_tar

echo "$input_files" \
    | xargs -I FILE bash -c 'retrieve_input "$@"' _ FILE

# Convert the FASTQ to a FASTA.
awk '
    NR % 4 == 1 { printf(">%s\n",substr($0,2)) }
    NR % 4 == 2 { print }
' "$output_fastq" > "$output_fasta"
