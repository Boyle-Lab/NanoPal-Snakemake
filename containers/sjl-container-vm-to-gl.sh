#!/usr/bin/env bash

set -euo pipefail

name="$1"

scp "vm:src/NanoPal-Snakemake/containers/${name}.sif" "${name}.sif"
scp "${name}.sif" "gl:/scratch/apboyle_root/apboyle0/slosh/nanopal-runs/containers/${name}.sif"
