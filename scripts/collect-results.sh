#!/usr/bin/env bash

set -euo pipefail

out_dir="$1"
shift

out_csv="$1"
shift

# Remaining arguments are the input directories.

# Pull the header line from an arbitrary CSV.
head -n 1 "$1"/result.csv > "$out_csv"

# Copy all the remaining CSV bodies, as well as the BED files of calls.
function copy_results {
    result_dir="$1"

    tail -n +2 "$result_dir"/result.csv >> "$out_csv"

    cp -t "$out_dir" "$result_dir"/nanopal_calls-*.bed
}

for result_dir in "$@"; do
    copy_results "$result_dir"
done





