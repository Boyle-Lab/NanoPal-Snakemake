#!/usr/bin/env bash

# Usage: ./run-glakes.sh runs/all.json …

set -euo pipefail

RUN="$1"
shift

export BLAST_USAGE_REPORT=0

snakemake \
  --use-singularity \
  --configfile "config/glakes.json" "$RUN" \
  --profile "profiles/glakes" \
  --snakefile "nanopal.smk" \
  "$@"


