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
        runtime="1h",
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
            "    --input {input.bam}"
            "    --workdir {params.workdir}"
            "    --output {wildcards.chromosome}"
            "    --ref_fa {input.reference}"
            "    --ref_ver {params.reference_version}"
            "    --type {params.mobile_element}"
            "    --chr {wildcards.chromosome}"
            "    --mode raw"
        )

rule gather_matches:
    # From Nanopal script:
    #
    # > pull out all reads having putative L1Hs signal reported by PALMER
    localrule: True
    log:
        scratch("_logs/gather_matches/{sample}.log"),
    benchmark:
        scratch("_benchmarks/gather_matches/{sample}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        palmer_sentinels=expand(
            scratch("palmer/{{sample}}/{chromosome}/done"),
            chromosome=config["chromosomes"],
        ),
        bam=scratch("alignment/{sample}/alignment.bam"),
    output:
        blast_matches=scratch("gather_matches/{sample}/blastn_refine.all.txt"),
        cigar_matches=scratch("gather_matches/{sample}/mapped.info.txt"),
    params:
        palmer_dir=scratch("palmer/{sample}/"),
    threads: 1
    shell:
        logged(
            "find {params.palmer_dir} -name 'blastn_refine.txt'"
            " | xargs cat"
            " > {output.blast_matches}"
            ,
            "samtools view {input.bam}"
            "   --min-MQ 10"
            "   --exclude-flags=0x100" # exclude secondary alignments
            "   --exclude-flags=0x200" # exclude not passed filters
            "   --exclude-flags=0x400" # exclude PCR/optical duplicates
            "   --exclude-flags=0x700" # exclude supplementary alignments
            " | awk '{{print $1, $3, $4, $6}}'" # read id, chromosome, position, cigar
            " > {output.cigar_matches}"
        )

rule parse_cigar:
    localrule: True
    log:
        scratch("_logs/parse_cigar/{sample}.log"),
    benchmark:
        scratch("_benchmarks/parse_cigar/{sample}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        container=containers("nanopal-binaries"),
        cigar_matches=scratch("gather_matches/{sample}/mapped.info.txt"),
    output:
        cigar_results=scratch("parse_cigar/{sample}/cigar_results.all.txt"),
        cigar_ref=scratch("parse_cigar/{sample}/cigar_ref.txt"),
        mapped_info=scratch("parse_cigar/{sample}/mapped.info.final.txt"),
    threads: 1
    shell:
        logged(
            "awk '{{print $4}}' {input.cigar_matches} | cigar_parser > {output.cigar_results}",
            "awk '{{print $4+$6+$10}}' {output.cigar_results} > {output.cigar_ref}",
            "paste {input.cigar_matches} {output.cigar_ref}"
            " | awk '{{print $1, $2, $3, $3+$5}}'"
            " > {output.mapped_info}"
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

rule _cigar:
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
        expand(scratch("parse_cigar/{sample}/mapped.info.final.txt"), sample=config["samples"])
