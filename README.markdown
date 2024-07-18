# Snakemake Pipeline for Nanopal

This is a port of <https://github.com/WeichenZhou/NanoPal-and-Cas9-targeted-enrichment-pipelines/>
from multiple bash scripts into a single Snakemake-orchestrated pipeline.

## Dependencies

To run, you'll need Snakemake and Singularity installed on the workstation or
cluster you want to use.

For Great Lakes you should be able to do `module load snakemake singularity` and
be all set.

For workstations, you'll need one with Singularity installed (`clover`
definitely has it, not sure about the others).  Then you'll need to install
Snakemake.  You could do this through Conda, or with a boring old virtualenv if
you prefer.

## Configuration

There are three layers of configuration.

`config/base.json` contains basic configuration for the pipeline itself, e.g.
which versions of the various tools to use.  There are also `config/local.json` and
`config/glakes.json` that let you change those values for locally running or
running on Great Lakes, but you likely won't need to do that (we should
probably just remove them entirely).

`profiles/glakes/config.yaml` contains the necessary configuration for using
Slurm on Great Lakes to parallelize the pipeline.  You might need to edit this
if you want to e.g. change which Slurm account the pipeline runs its tasks
under.

`runs/*.json` is where you'll configure an individual run of a pipeline.
A single run config here looks something like this:

```json
{
    "batch_id": "real-data-test-02",
    "mobile_elements": ["LINE", "AluYa", "AluYb"],
    "reference_version": "GRCh38",
    "reference":      "/scratch/apboyle_root/apboyle0/slosh/nanopal-snakemake-test-runs/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    "data_path":      "/nfs/turbo/boylelab/nanopore_data/MEI/AD",
    "scratch_path":   "/scratch/apboyle_root/apboyle0/slosh/nanopal-snakemake-test-runs/scratch",
    "container_path": "/scratch/apboyle_root/apboyle0/slosh/nanopal-snakemake-test-runs/containers",
    "samples": [
        "20230202_1513_MN35288_FAS89555_2289ca4f"
    ],
    "chromosomes": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ]
}
```

`batch_id` is a name for the run.  The pipeline will store all its output under
a directory with this name in the `scratch_path` directory.

`data_path` is the path that contains the directories of your samples.  Each of
`samples` should be a directory inside that data path, which contains
a `basecalled_output.tar` file.  Note that some of our older sequencing runs
have a slightly different format here — see the TODO in
`scripts/retrieve-input.sh` for more context.

`container_path` is the path to a directory to store the Singularity containers
(see below).

## Containers

The pipeline will download its own containers to the directory specified in
`container_path`, *except* for the custom containers that we need to build
ourselves.  For those you'll need to build them yourself (or find someone who's
got them built) and copy them over manually before you run for the first time.

Unfortunately you can't build Singularity containers on Great Lakes, so you'll
need to have Singularity installed on a machine you have `root` access to (or
that's been specially configured to allow non-`root` container building). Once
you've got that you should be able to build the `*.sif` files with:

    singularity build --fakeroot nanopal-binaries.sif nanopal-binaries.def
    singularity build --fakeroot palmer.sif           palmer.def

Then copy those to the `container_path` wherever you're trying to run.

There may be some other snags I've forgotten to document here, ping me (`slosh`)
if you get an error here.

## Running

Once you've got your `runs/whatever.json` and the manually-built containers
ready, you should be all set to run.  There are a couple of helper scripts that
will make it a little less tedious to invoke Snakemake:

* `run-local.sh` for running directly on a workstation.
* `run-glakes.sh` for running on Great Lakes via Slurm.

Take a look at the individual scripts to see exactly what they're doing, but
here are some examples to get you started.  To do a dry run and just print what
Snakemake thinks it needs to do:

    script          run config                        snakemake args
    vvvvvv          vvvvvvvvvv                        vvvvvvvvvvvvvv
    ./run-glakes.sh samples/real-data-test-gl.json    _all --cores=36 --dry-run
    ./run-local.sh  samples/real-data-test-local.json _all --cores=12 --dry-run

Unfortunately you *do* need to pass `--cores` even for a dry run.  The reason is
that some of the commands that get invoked depend in a non-trivial way on the
number of cores used, which means Snakemake will consider them changed if the
number of cores changes (and the default `--cores` is `1`), which will make the
dry run look like it'll have to do more work that the real run actually will.

To do a full run for real on a workstation (e.g. `clover`):

    ./run-local.sh samples/real-data-test-local.json _all --cores=12 > snakemake.log

Note that you'll probably want to do this in a `screen` or `tmux` session so it
doesn't get killed if your wifi dies.

When running on the cluster you're not supposed to do any real work on the login
nodes, so you should use Slurm to invoke the `./run-glakes.sh` call on a Slurm
worker node.  The `run.sbat` file does that, so copy and edit that to point at
your run config and use `sbatch my-run.sbat` to kick off the Snakemake
controller, which will then fan out the tasks to other nodes.

## Output

The pipeline will put all of its output into subdirectories of `<scratch_path>/<batch_id>`:

```
/scratch/apboyle_root/apboyle0/slosh/nanopal-snakemake-test-runs/scratch
└── [ 4.0K]  real-data-test-02/
    ├── [ 4.0K]  _benchmarks/
    ├── [ 4.0K]  _logs/
    ├── [ 4.0K]  alignment/
    ├── [ 4.0K]  find_on_target/
    ├── [ 4.0K]  find_valid_reads/
    ├── [ 4.0K]  gather_matches/
    ├── [ 4.0K]  input/
    ├── [ 4.0K]  intersect/
    ├── [ 4.0K]  intersect_again/
    ├── [ 4.0K]  palmer/
    ├── [ 4.0K]  palmer_on_target/
    ├── [ 4.0K]  parse_cigar/
    └── [ 4.0K]  reference/
```

The `_benchmarks` directory stores the benchmark `tsv`s collected automatically
by Snakemake, which can be nice if you want to try to optimize the
`runtime`/`memory` requested for a particular task.

`_logs` stores the logs for all the various Snakemake tasks.

All the rest of the directories are the working directories for the Snakemake
tasks with the various names.  Each will have one subdirectory per sample, and
those that are parallelized across targets will have another layer of
subdirectories under them:

```
find_valid_reads/
└── [ 4.0K]  20230202_1513_MN35288_FAS89555_2289ca4f/
    └── [ 5.5M]  RC.all.list

parse_cigar/
└── [ 4.0K]  20230202_1513_MN35288_FAS89555_2289ca4f/
    ├── [ 4.0K]  AluYa/
    │   ├── [ 1.3M]  cigar_ref.txt
    │   ├── [  13M]  cigar_results.all.txt
    │   └── [  18M]  mapped.info.final.txt
    ├── [ 4.0K]  AluYb/
    │   ├── [ 1.3M]  cigar_ref.txt
    │   ├── [  13M]  cigar_results.all.txt
    │   └── [  18M]  mapped.info.final.txt
    └── [ 4.0K]  LINE/
        ├── [ 1.3M]  cigar_ref.txt
        ├── [  13M]  cigar_results.all.txt
        └── [  18M]  mapped.info.final.txt
```

The final output numbers (which were just dumped to stdout in the vanilla
version) will be under `intersect_again/result-log.txt` (TODO we should add
another Snakemake step to collect all the end results into a single place for
easier viewing/downloading).

TODO Note, at some point we might invert the structure of the directories so
instead of `<step>/<sample>/<target>` it's structured
`<sample>/<step>/<target>`, but that's still yet to be decided.
