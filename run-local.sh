#!/usr/bin/env bash

# Usage: ./run-local.sh runs/all.json …

set -euo pipefail

RUN="$1"
shift

export BLAST_USAGE_REPORT=0

snakemake \
  --use-singularity \
  --configfile "config/local.json" "$RUN" \
  --snakefile "nanopal.smk" \
  --singularity-args "--bind /home/slosh/scratch/" \
  "$@"
