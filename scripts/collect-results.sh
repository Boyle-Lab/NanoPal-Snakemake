#!/usr/bin/env bash

set -euo pipefail

out_dir="$1"
shift

# Remaining arguments are the input directories.

out_csv="${out_dir}/results.csv"


# Pull the header line from an arbitrary CSV.
head -n 1 "$1"/result.csv > "$out_csv"


# Copy all the remaining CSV bodies.
function copy_results {
    tail -n +2 "${1}/result.csv" >> "$out_csv"
}

for result_dir in "$@"; do
    copy_results "$result_dir"
done





