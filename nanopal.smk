import os


# Config ---------------------------------------------------------------------

configfile: "config/base.json"

SAMPLES = config["samples"]
BATCH_ID = config["batch_id"]
DATA_PATH = config["data_path"]
SCRATCH_PATH = config["scratch_path"]
CONTAINER_PATH = config["container_path"]

wildcard_constraints:
    id="[^/]+",


# Utilities -------------------------------------------------------------------

def data(path):
    return os.path.join(DATA_PATH, path)

def scratch(path):
    return os.path.join(SCRATCH_PATH, BATCH_ID, path)

def containers(name):
    return os.path.join(CONTAINER_PATH, f"{name}.sif")

def logged(shell_command):
    return f"time ({shell_command}) " + ">{log:q} 2>&1"


# Rules -----------------------------------------------------------------------

rule container:
    localrule: True
    output:
        containers('{name}')
    params:
        source=lambda wc: config["containers"][wc.name],
    shell:
        "singularity pull {output} {params.source}"



# PHONY -----------------------------------------------------------------------

rule _containers:
    localrule: True
    input:
        expand(containers('{name}'), name=config["containers"].keys())
