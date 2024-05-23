#!/usr/bin/env bash

# Usage: ./run-glakes-slurm.sh samples/all.json …

set -euo pipefail

SAMPLES="$1"
shift

snakemake \
  --use-singularity \
  --configfile "config/glakes-slurm.json" "$SAMPLES" \
  --profile "profiles/glakes-slurm" \
  "$@"
