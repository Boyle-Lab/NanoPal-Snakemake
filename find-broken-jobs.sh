#!/usr/bin/env bash

set -euo pipefail

slurm_log="$1"

awk '
BEGIN {
    name = snakemake_id = 0
    OFS=","
}

/^rule .*:$/ {
    match($0, /rule ([^ :]+)/, m)
    name = m[1]
}

/ +jobid: [0-9]+/ {
    snakemake_id = $2
}

/^$/ {
    if (name) {
        job_names[snakemake_id] = name
        name = snakemake_id = 0
    }
}

/Submitted job [0-9]+ with external jobid/ {
    match($0, /Submitted job ([0-9]+) with external jobid.*batch job ([0-9]+)/, m)
    snakemake_id = m[1]
    slurm_id = m[2]
    slurm_ids[snakemake_id] = slurm_id
}

/^Finished job/ {
    match($0, /^Finished job ([0-9]+)/, m)
    finished_jobs[m[1]] = 1
}

END {
    print "job_name", "snakemake_id", "slurm_id", "log_file_pattern"

    for (snakemake_id in slurm_ids) {
        if (!finished_jobs[snakemake_id]) {
            name = job_names[snakemake_id]
            slurm_id = slurm_ids[snakemake_id]
            log_file = "slurm-logs/" name "/" name "-*-" slurm_id ".out"
            print name, snakemake_id, slurm_id, log_file
        }
    }
}
' "$slurm_log" | column -s, -t

# e.g.:
#
# ./find-broken-jobs.sh snakemake-batch-run-49977560.out | f 4 | tail -n +2 | xargs -n1 sh -c 'cat $1; echo XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX' -- | less
