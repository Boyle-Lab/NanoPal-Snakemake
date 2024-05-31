#!/usr/bin/env bash

set -eu

palmer_reads="$1" # read.all.palmer.final.txt
out_dir="$2"

mkdir -p "$out_dir"
cd "$out_dir"

awk '{print $11, $12, $12+1, $1}' "$palmer_reads"  | grep -v "NON" > read.5.Q.txt
awk '{print $11, $13, $13+1, $1}' "$palmer_reads"  | grep -v "NON" > read.3.Q.txt

cp read.5.Q.txt Q.txt
