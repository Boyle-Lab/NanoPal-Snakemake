#!/usr/bin/env bash

set -euo pipefail

out_csv="$1"
shift

# Remaining arguments are the input CSVs.

# Pull the header line from an arbitrary CSV.
head -n 1 "$1" > "$out_csv"


# Copy all the remaining CSV bodies.
function copy_results {
    tail -n +2 "$1" >> "$out_csv"
}

for result_csv in "$@"; do
    copy_results "$result_csv"
done





