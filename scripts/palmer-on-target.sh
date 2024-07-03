#!/usr/bin/env bash

set -eu

# Input paths
palmer_blast="$1" # blastn_refine.all.txt
nanopal_reads="$2" # read.all.txt
palmer_map="$3" # mapped.info.final.txt
mei_cut_site="$4"

function full_input() {
    # Join the PALMER BLAST results with the Nanopal blast, then join with the
    # alignment mappings.  Output as:
    #
    #     read_id read_length
    #     nanopal_plus5 nanopal_minus5 nanopal_plus3 nanopal_minus3
    #     qacc evalue qstart qend sstart send
    #     chromosome align_start_pos align_end_pos
    join -1 1 -2 2 \
         -a 1 -a 2 \
         -o auto -e 'null' \
         --check-order \
         <(cat "$nanopal_reads" | sort -k 1) \
         <(cat "$palmer_blast" | sed -Ee 's/_[_0-9]*//' | sort -k 2) \
    | join -1 1 -2 1 \
           -a 1 -a 2 \
           -o auto -e 'null' \
           --check-order \
           - \
           <(cat "$palmer_map" | sort -k 1)
}

awk \
    -v cut_lo=100 \
    -v cut_hi="$mei_cut_site" \
    '
    BEGIN { OFS = "\t" }

    function print_result() {
        if (chr == "null") { chr = "NON" }
        if (start_pos == "null") { start_pos = 0 }
        if (end_pos == "null") { end_pos = 0 }

        print read_name, read_length, n_plus5, n_minus5, n_plus3, n_minus3, p_plus5, p_minus5, p_plus3, p_minus3, chr " " start_pos " " end_pos
    }

    {
        if (read_name != $1) {
            if (read_name != "") {
                print_result()
            }

            read_name = $1
            read_length = $2
            len_fix = read_length - 100

            p_plus5 = 0
            p_minus5 = 0
            p_plus3 = 0
            p_minus3 = 0
        }

        # These should all be the same for all the rows, so overwriting them
        # every time should be fine.
        n_plus5   = $3
        n_minus5  = $4
        n_plus3   = $5
        n_minus3  = $6

        # These are the things changing for a single read (since it might have
        # multiple BLAST results).
        qacc      = $7
        evalue    = $8
        qstart    = $9
        qend      = $10
        sstart    = $11
        send      = $12

        # These should also all be the same for all the rows.
        chr       = $13
        start_pos = $14
        end_pos   = $15

        if (evalue != "null") {
            # These BLAST results come from PALMER, which reverses query/subject
            # compared to the other script here, so we do not need to do the
            # extra swap the other script did.

            if(qend > cut_hi) {
                if (sstart < cut_lo || send < cut_lo) {
                    if (sstart < send) {
                        p_plus5 = 1
                    } else {
                        p_minus5 = 1
                    }
                }

                if (sstart > len_fix || send > len_fix) {
                    if (sstart < send) {
                        p_plus3 = 1
                    } else {
                        p_minus3 = 1
                    }
                }
            }
        }
    }

    END { if (read_name) { print_result() } }
' <(full_input)
