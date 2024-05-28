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


def logged(*shell_commands):
    commands = "\n".join(shell_commands)
    return f"time ({commands}) " + ">{log:q} 2>&1"


# Rules -----------------------------------------------------------------------


rule container:
    localrule: True
    output:
        containers("{name}"),
    params:
        source=lambda wc: config["containers"][wc.name],
    shell:
        "singularity pull {output} {params.source}"


rule input:
    localrule: True
    log:
        scratch("_logs/input/{sample}.log"),
    input:
        data("{sample}"),
    output:
        fastq=scratch("input/{sample}/batch.fastq"),
        fasta=scratch("input/{sample}/batch.fasta"),
    threads: 1
    shell:
        logged(
            "cat {input}/*.fastq > {output.fastq}",
            r"""cat {output.fastq} | awk '
                                        NR % 4 == 1 {{ printf(">%s\n",substr($0,2)) }}
                                        NR % 4 == 2 {{ print }}
                                    ' > {output.fasta}""",
        )


rule index_reference:
    log:
        scratch("_logs/index_reference.log"),
    benchmark:
        scratch("_benchmarks/index_reference.tsv")
    container:
        containers("minimap2")
    input:
        container=containers("minimap2"),
        reference=config["reference"],
    output:
        scratch("reference/index.mmi"),
    threads: 16
    resources:
        mem="24GB",
        runtime="30m",
    shell:
        logged(
            "minimap2"
            " -a"
            " -x map-ont"
            " -d {output}"
            " -t {threads}"
            " {input.reference}"
        )


rule alignment:
    log:
        scratch("_logs/alignment/{sample}.log"),
    benchmark:
        scratch("_benchmarks/alignment/{sample}.tsv")
    container:
        containers("minimap2")
    input:
        container=containers("minimap2"),
        fastq=scratch("input/{sample}/batch.fastq"),
        index=scratch("reference/index.mmi"),
    output:
        scratch("alignment/{sample}/alignment.sam"),
    threads: 16
    resources:
        mem="48GB",
        runtime="30m",
    shell:
        logged(
            "minimap2"
            " -a"
            " -x map-ont"
            " -t {threads}"
            " {input.index} {input.fastq}"
            " > {output}"
        )


rule index_alignment:
    log:
        scratch("_logs/index_alignment/{sample}.log"),
    benchmark:
        scratch("_benchmarks/index_alignment/{sample}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        sam=scratch("alignment/{sample}/alignment.sam"),
    output:
        bam=scratch("alignment/{sample}/alignment.bam"),
        bai=scratch("alignment/{sample}/alignment.bam.bai"),
    threads: 2  # TODO
    resources:
        mem="8GB",
        runtime="20m",
    shell:
        logged(
            r"""
            samtools sort -o {output.bam} {input.sam} 
            samtools index {output.bam}
            """
        )


rule palmer:
    log:
        scratch("_logs/palmer/{sample}_{chromosome}.log"),
    benchmark:
        scratch("_benchmarks/palmer/{sample}_{chromosome}.tsv")
    container:
        containers("palmer")
    input:
        container=containers("palmer"),
        reference=config["reference"],
        bam=scratch("alignment/{sample}/alignment.bam"),
    output:
        sentinel=touch(scratch("palmer/{sample}/{chromosome}/done")),
    params:
        workdir=scratch("palmer/{sample}/{chromosome}/"),
        reference_version=config["reference_version"],
        mobile_element=config["mobile_element"],
    threads: 2  # TODO
    resources:
        mem="4GB",
        runtime="15m",
    shell:
        logged(
            "/palmer/PALMER"
            "  --input {input.bam}"
            "  --workdir {params.workdir}"
            "  --output {wildcards.chromosome}"
            "  --ref_fa {input.reference}"
            "  --ref_ver {params.reference_version}"
            "  --type {params.mobile_element}"
            "  --chr {wildcards.chromosome}"
            "  --mode raw"
        )


# PHONY -----------------------------------------------------------------------


rule _containers:
    localrule: True
    input:
        expand(containers("{name}"), name=config["containers"].keys()),


rule _input:
    localrule: True
    input:
        expand(scratch("input/{sample}/batch.fastq"), sample=config["samples"]),
        expand(scratch("input/{sample}/batch.fasta"), sample=config["samples"]),


rule _alignment:
    localrule: True
    input:
        expand(scratch("alignment/{sample}/alignment.bam"), sample=config["samples"]),


rule _palmer:
    localrule: True
    input:
        expand(
            scratch("palmer/{sample}/{chromosome}/done"),
            sample=config["samples"],
            chromosome=config["chromosomes"],
        ),


rule _all:
    localrule: True
    input:
        rules._palmer.input,
        rules._alignment.input,
