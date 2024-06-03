#!/usr/bin/env bash

set -eu

bam="$1" # Nanopore.sorted.bam
ref_mei=$(readlink -f "$2") # hg38.RM.L1.ref
in_summary="$3" # summary.final.txt

out_dir="$4"
out_summary="$5" # summary.final.2.txt

mkdir -p "$out_dir"
cd "$out_dir"

##different types of LINE
#######MEI
grep    L1PA "$ref_mei" > ref.L1PA
grep    L1HS "$ref_mei" > ref.L1HS
grep -v L1HS "$ref_mei"| grep -v L1PA | grep L1 > ref.L1

##Valid mapping reads
samtools view "$bam" -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 -f 0x10 \
    | awk '{print $1}' \
    >  RC.all.list

rm -f RC.tag
while read -r aa b c d e f g h i j  k l m n o p q r
do
awk -v aa="$aa" '
    $1 == aa { print 1; exit 0 }
    END { print 0 }
' RC.all.list >> RC.tag

wait
done < "$in_summary"

paste "$in_summary" RC.tag > "$out_summary"
