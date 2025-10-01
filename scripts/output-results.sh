#!/usr/bin/env bash

set -euo pipefail
set -x

dataset_id="$1"; shift
mei="$1"; shift # LINE
bam="$1"; shift # alignment.bam
potential_meis="$1"; shift # potential.clustered.txt.fi
result_log="$1"; shift # result-log.txt
threads="$1"; shift
out_dir="$1"; shift
out_csv="$1"; shift
out_log="$1"; shift
out_bed="$1"; shift
out_bed_multi="$1"; shift

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
potential_mei_count=$(awk 'BEGIN {n=0}        {n+=1} END {print n}' "$potential_meis")
multiple_support=$(awk    'BEGIN {n=0} $4 > 1 {n+=1} END {print n}' "$potential_meis")

# These statistics take a while to compute.  Unfortunately trying to parallelize
# it with --threads to samtools and/or seqkit doesn't really help all that much.
# But we can at least parallelize the individual stat computations.

# Get some basic read counts/stats from the BAM.
total_reads=$(samtools view --threads "$threads" -c "$bam") # need to compute first to divide by in the next stats
samtools view --threads 2 -c -F 0x4 "$bam" | awk -v n="$total_reads" '{print 100.0 * $1/n}'         > "$out_dir"/stat_percent_mapped_reads &
samtools view --threads 2 -c -f 0x4 "$bam" | awk -v n="$total_reads" '{print 100.0 * $1/n}'         > "$out_dir"/stat_percent_unmapped_reads &
samtools view --threads 2 -c -q 10  "$bam" | awk -v n="$total_reads" '{print 100.0 * (1 - ($1/n))}' > "$out_dir"/stat_percent_low_mapq &

# Get some stats about the unmapped reads. seqkit stats output looks like:
#
#     file    format  type    num_seqs        sum_len min_len avg_len max_len Q1      Q2      Q3      sum_gap N50     Q20(%)  Q30(%)  GC(%)
#     data…

samtools fastq -F 0x4 "$bam" --threads 1 | seqkit stats -a -T --threads 1 | awk 'NR == 2 {print  $7}' > "$out_dir"/stat_mapped_average_length &
samtools fastq -F 0x4 "$bam" --threads 1 | seqkit stats -a -T --threads 1 | awk 'NR == 2 {print $13}' > "$out_dir"/stat_mapped_n50 &
samtools fastq -f 0x4 "$bam" --threads 1 | seqkit stats -a -T --threads 1 | awk 'NR == 2 {print  $7}' > "$out_dir"/stat_unmapped_average_length &
samtools fastq -f 0x4 "$bam" --threads 1 | seqkit stats -a -T --threads 1 | awk 'NR == 2 {print $13}' > "$out_dir"/stat_unmapped_n50 &

wait

percent_mapped_reads=$(cat   "$out_dir"/stat_percent_mapped_reads)
percent_unmapped_reads=$(cat "$out_dir"/stat_percent_unmapped_reads)
percent_low_mapq=$(cat       "$out_dir"/stat_percent_low_mapq)

mapped_average_length=$(cat   "$out_dir"/stat_mapped_average_length)
mapped_n50=$(cat              "$out_dir"/stat_mapped_n50)
unmapped_average_length=$(cat "$out_dir"/stat_unmapped_average_length)
unmapped_n50=$(cat            "$out_dir"/stat_unmapped_n50)

# Output CSV ------------------------------------------------------------------

{
    echo "id,target,total_reads,percent_mapped_reads,mapped_average_length,mapped_n50,percent_unmapped_reads,unmapped_average_length,unmapped_n50,percent_low_mapq,signal_single_end,signal_double_end,no_signal,enrichment,potential_meis,multiply_supported_potential_meis"

    echo "${dataset_id},${mei},${total_reads},${percent_mapped_reads}%,${mapped_average_length},${mapped_n50},${percent_unmapped_reads}%,${unmapped_average_length},${unmapped_n50},${percent_low_mapq}%,${signal1},${signal2},${fail},${enrichment},${potential_mei_count},${multiple_support}"
} > "$out_csv"

cp "$result_log" "$out_log"

# Output BEDs -----------------------------------------------------------------

awk '
  BEGIN {
    FS="\t"
    OFS="\t"
  }
  {
    print $1, $2, $3, $4 " " $6
  }
' "$potential_meis" > "$out_bed"

awk '$4 > 1' "$out_bed" > "$out_bed_multi"
