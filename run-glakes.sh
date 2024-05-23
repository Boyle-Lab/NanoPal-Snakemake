#!/usr/bin/env bash

# Usage: ./run-glakes.sh samples/all.json …

set -euo pipefail

SAMPLES="$1"
shift

snakemake \
  --use-singularity \
  --configfile "config/glakes.json" "$SAMPLES" \
  --profile "profiles/glakes" \
  --snakefile "nanopal.smk"
  "$@"
