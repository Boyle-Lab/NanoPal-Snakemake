#!/usr/bin/env bash

# Usage: ./run-local.sh samples/all.json …

set -euo pipefail

SAMPLES="$1"
shift

export BLAST_USAGE_REPORT=0

snakemake \
  --use-singularity \
  --configfile "config/local.json" "$SAMPLES" \
  --snakefile "nanopal.smk" \
  --singularity-args "--bind /home/slosh/scratch/" \
  "$@"
