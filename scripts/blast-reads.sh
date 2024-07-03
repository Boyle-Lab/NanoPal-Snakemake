#!/usr/bin/env bash

set -euo pipefail

mei_db="$1"
reads_fasta="$2"
mei_cut_site="$3"
blast_threads="$4"

function read_lengths {
    # Compute the length of each read in the FASTA, output as:
    #
    #     id length
    awk '
        /^>/ {
            if (name) {
                print substr(name, 2), n
            }
            name = $1
            n = 0
        }
        /^[^>]/ {
            n += length($0)
        }
        END {
            if (name) {
                print substr(name, 2), n
            }
        }
    ' "$reads_fasta" \
        | sort -k 1
}

function blast_query {
    # BLAST the query (i.e. the reads) against the consensus sequence.  Output
    # as (qacc = query accession = read id):
    #
    #     qacc sacc evalue qstart qend sstart send

    blastn \
        -evalue 0.001 \
        -task blastn \
        -gapopen 4 \
        -num_threads "$blast_threads" \
        -query "$reads_fasta" \
        -db "$mei_db" \
        -outfmt "6 qacc sacc evalue qstart qend sstart send" \
        | sort -k 1
}

function full_input {
    # Join the BLAST results with the read lengths.  Output as:
    #
    #     id length qacc evalue qstart qend sstart send
    join -1 1 -2 1 -a 1 -a 2 -o auto -e 'null' \
        --check-order \
        <(read_lengths) <(blast_query)
}


# Compute hits.
awk \
    -v cut_lo=100 \
    -v cut_hi="$mei_cut_site" \
    '
    function print_result() {
        printf("%s\t%s\t%d\t%d\t%d\t%d\n",
               read_name, read_length, plus5, minus5, plus3, minus3)
    }

    {
        if (read_name != $1) {
            if (read_name != "") {
                print_result()
            }

            read_name = $1
            read_length = $2
            len_fix = read_length - 100
            plus5 = 0
            minus5 = 0
            plus3 = 0
            minus3 = 0
        }

        qacc   = $3
        evalue = $4
        qstart = $5
        qend   = $6
        sstart = $7
        send   = $8

        if (evalue != "null") {
            # BLAST reports results from the "perspective" of the query, i.e.
            # you will always have qstart < qend, and then sstart and send will
            # depend on the direction of the match.  This is the opposite of
            # what we want here because we are using the L1 as the subject, and
            # we want to look for the cut site location in THAT first.  To deal
            # with this we can detect when we have a swap and flip both sets of
            # positions, to make the READS the ones that are direction-dependent
            # again.
            if (sstart > send) {
                tmp = qend; qend = qstart; qstart = tmp
                tmp = send; send = sstart; sstart = tmp
            }

            if(send > cut_hi) {
                if (qstart < cut_lo || qend < cut_lo) {
                    if (qstart < qend) {
                        plus5 = 1
                    } else {
                        minus5 = 1
                    }
                }

                if (qstart > len_fix || qend > len_fix) {
                    if (qstart < qend) {
                        plus3 = 1
                    } else {
                        minus3 = 1
                    }
                }
            }
        }
    }

    END { if (read_name) { print_result() } }
' <(full_input)

