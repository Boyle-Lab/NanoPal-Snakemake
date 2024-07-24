#!/usr/bin/env bash

set -euo pipefail
set -x

dataset_id="$1"
mei="$2" # LINE
bam="$3" # alignment.bam
potential_meis="$4" # potential.clustered.txt.fi
result_log="$5" # result-log.txt
out_csv="$6"
out_log="$7"

# Compute statistics ----------------------------------------------------------

# This is a little ugly, sorry, but we'll parse the result log for these:
#
#     There are ${SIGNAL1} reads capturing putative signals on one end.
#     There are ${SIGNAL2} reads capturing putative signals on both ends.
#     There are ${FAIL} reads having no putative signals.
#     Enrichment: ${ENRICHMENT}% of reads have signal(s).

signal1=$(awk    'NR == 1 {print $3}' "$result_log")
signal2=$(awk    'NR == 2 {print $3}' "$result_log")
fail=$(awk       'NR == 3 {print $3}' "$result_log")
enrichment=$(awk 'NR == 4 {print $2}' "$result_log")

# Count potential events with/without multiple (≥2) read support.
potential_mei_count=$(awk '          {n+=1} END {print n}' "$potential_meis")
multiple_support=$(awk    '$5+$6 > 1 {n+=1} END {print n}' "$potential_meis")

# Get some basic read counts/stats from the BAM.
total_reads=$(samtools    view -c        "$bam")
percent_mapped_reads=$(samtools   view -c -F 0x4 "$bam" | awk -v n="$total_reads" '{print $1/n}')
percent_unmapped_reads=$(samtools view -c -f 0x4 "$bam" | awk -v n="$total_reads" '{print $1/n}')
percent_low_mapq=$(samtools view -c -q 10 "$bam" | awk -v n="$total_reads" '{print 1 - ($1/n)}')

# Get some stats about the unmapped reads. seqkit stats output looks like:
#
#     file    format  type    num_seqs        sum_len min_len avg_len max_len Q1      Q2      Q3      sum_gap N50     Q20(%)  Q30(%)  GC(%)
#     data…

unmapped_average_length=$(samtools view -f 0x4 "$bam" | samtools fastq | seqkit stats -a -T | awk 'NR == 2 {print  $7}')
unmapped_n50=$(samtools            view -f 0x4 "$bam" | samtools fastq | seqkit stats -a -T | awk 'NR == 2 {print $13}')

# Output CSV ------------------------------------------------------------------

{
    echo "id,target,total_reads,percent_mapped_reads,percent_unmapped_reads,unmapped_average_length,unmapped_n50,percent_low_mapq,signal_single_end,signal_double_end,no_signal,enrichment,potential_meis,multiply_supported_potential_meis"

    echo "${dataset_id},${mei},${total_reads},${percent_mapped_reads},${percent_unmapped_reads},${unmapped_average_length},${unmapped_n50},${percent_low_mapq},${signal1},${signal2},${fail},${enrichment},${potential_mei_count},${multiple_support}"
} > "$out_csv"

cp "$result_log" "$out_log"
