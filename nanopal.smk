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
        script="scripts/retrieve-input.sh",
        input_dir=data("{sample}"),
    output:
        fastq=scratch("input/{sample}/batch.fastq"),
        fasta=scratch("input/{sample}/batch.fasta"),
    threads: 1
    shell:
        logged(
            "./{input.script} "
            "  {input.input_dir}"
            "  {output.fastq}"
            "  {output.fasta}"
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
    threads: 2
    resources:
        mem="16GB",
        runtime="10m",
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
        mem="32GB",
        runtime="30m",
    shell:
        logged(
            "minimap2"
            " -a"
            " -x map-ont"
            " -t {threads}"
            " --secondary=no"
            " --eqx"
            " -Y"
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

rule find_valid_reads:
    localrule: True
    log:
        scratch("_logs/find_valid_reads/{sample}.log"),
    benchmark:
        scratch("_benchmarks/find_valid_reads/{sample}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        bam=scratch("alignment/{sample}/alignment.bam"),
    output:
        valid_read_ids=scratch("find_valid_reads/{sample}/RC.all.list"),
    threads: 2
    shell:
        logged(
            "samtools view {input.bam}"
            "    -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 -f 0x10"
            " | awk '{{print $1}}'"
            " > {output.valid_read_ids}"
        )

def palmer_mei_param(wc):
    return {
        'LINE': 'LINE',
        'AluYa': 'ALU',
        'AluYb': 'ALU',
    }[wc.mei]

rule palmer:
    log:
        scratch("_logs/palmer/{sample}_{mei}_{chromosome}.log"),
    benchmark:
        scratch("_benchmarks/palmer/{sample}_{mei}_{chromosome}.tsv")
    container:
        containers("palmer")
    input:
        container=containers("palmer"),
        reference=config["reference"],
        bam=scratch("alignment/{sample}/alignment.bam"),
    output:
        sentinel=touch(scratch("palmer/{sample}/{mei}/{chromosome}/done")),
    params:
        workdir=scratch("palmer/{sample}/{mei}/{chromosome}/"),
        reference_version=config["reference_version"],
        palmer_mei=palmer_mei_param,
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
            "    --type {params.palmer_mei}"
            "    --chr {wildcards.chromosome}"
            "    --mode raw"
        )

rule gather_matches:
    # From Nanopal script:
    #
    # > pull out all reads having putative L1Hs signal reported by PALMER
    localrule: True
    log:
        scratch("_logs/gather_matches/{sample}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/gather_matches/{sample}_{mei}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        script="scripts/gather-palmer-results.sh",
        palmer_sentinels=expand(
            scratch("palmer/{{sample}}/{{mei}}/{chromosome}/done"),
            chromosome=config["chromosomes"],
        ),
        bam=scratch("alignment/{sample}/alignment.bam"),
    output:
        palmer_blast=scratch("gather_matches/{sample}/{mei}/blastn_refine.all.txt"),
        palmer_cigar=scratch("gather_matches/{sample}/{mei}/mapped.info.txt"),
    params:
        palmer_dir=scratch("palmer/{sample}/{mei}/"),
    threads: 1
    shell:
        logged(
            "./{input.script} "
            "  {params.palmer_dir}"
            "  {input.bam}"
            "  {output.palmer_blast}"
            "  {output.palmer_cigar}"
        )

rule parse_cigar:
    localrule: True
    log:
        scratch("_logs/parse_cigar/{sample}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/parse_cigar/{sample}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        container=containers("nanopal-binaries"),
        cigar_matches=scratch("gather_matches/{sample}/{mei}/mapped.info.txt"),
    output:
        cigar_results=scratch("parse_cigar/{sample}/{mei}/cigar_results.all.txt"),
        cigar_ref=scratch("parse_cigar/{sample}/{mei}/cigar_ref.txt"),
        mapped_info=scratch("parse_cigar/{sample}/{mei}/mapped.info.final.txt"),
    threads: 1
    shell:
        logged(
            "awk '{{print $4}}' {input.cigar_matches} | cigar_parser > {output.cigar_results}",
            "awk '{{print $4+$6+$10}}' {output.cigar_results} > {output.cigar_ref}",
            "paste {input.cigar_matches} {output.cigar_ref}"
            " | awk '{{print $1, $2, $3, $3+$5}}'"
            " > {output.mapped_info}"
        )

def mei_fasta(wc):
    return {
        'LINE': "meis/L1.3",
        'AluYa': "meis/AluYa5",
        'AluYb': "meis/AluYb8",
    }[wc.mei]

def mei_cut_site(wc):
    return {
        'LINE': 5900,
        'AluYa': 225,
        'AluYb': 227,
    }[wc.mei]

rule find_on_target:
    log:
        scratch("_logs/find_on_target/{sample}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/find_on_target/{sample}_{mei}.tsv")
    container:
        containers("blast")
    input:
        container=containers("blast"),
        script="scripts/blast-reads.sh",
        reads_fasta=scratch("input/{sample}/batch.fasta"),
        mei_fasta=mei_fasta,
    output:
        nanopal_reads=scratch("find_on_target/{sample}/{mei}/read.all.txt")
    params:
        mei_cut_site=mei_cut_site,
    threads: 1
    resources:
        mem="16GB",
        runtime="15m",
    shell:
        logged(
            "./{input.script}"
            "  {input.mei_fasta}"
            "  {input.reads_fasta}"
            "  {params.mei_cut_site}"
            "  > {output}"
        )

rule palmer_on_target:
    localrule: True
    log:
        scratch("_logs/palmer_on_target/{sample}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/palmer_on_target/{sample}_{mei}.tsv")
    container:
        containers("palmer")
    input:
        container=containers("palmer"),
        script="scripts/palmer-on-target.sh",
        palmer_blast=scratch("gather_matches/{sample}/{mei}/blastn_refine.all.txt"),
        palmer_map=scratch("parse_cigar/{sample}/{mei}/mapped.info.final.txt"),
        nanopal_reads=scratch("find_on_target/{sample}/{mei}/read.all.txt"),
    output:
        scratch("palmer_on_target/{sample}/{mei}/read.all.palmer.final.txt")
    params:
        mei_cut_site=mei_cut_site,
    threads: 2
    shell:
        logged(
            "./{input.script}"
            "  {input.palmer_blast}"
            "  {input.nanopal_reads}"
            "  {input.palmer_map}"
            "  {params.mei_cut_site}"
            "  > {output}"
        )

def ref_mei(wc):
    return {
        "LINE":  "meis/hg38.RM.L1.ref",
        "AluYa": "meis/hg38.RM.ALU.ref",
        "AluYb": "meis/hg38.RM.ALU.ref",
    }[wc.mei]

def orig_mei(wc):
    return {
        "LINE":  "meis/PALMER.NA12878.L1.txt",
        "AluYa": "meis/PALMER.NA12878.ALU.txt",
        "AluYb": "meis/PALMER.NA12878.ALU.txt",
    }[wc.mei]

rule intersect:
    log:
        scratch("_logs/intersect/{sample}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/intersect/{sample}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        script="scripts/intersect.sh",
        container=containers("nanopal-binaries"),
        palmer_reads=scratch("palmer_on_target/{sample}/{mei}/read.all.palmer.final.txt"),
        ref_mei=ref_mei,
        orig_mei=orig_mei,
    output:
        out_dir=directory(scratch("intersect/{sample}/{mei}/")),
        out_summary=scratch("intersect/{sample}/{mei}/summary.final.txt"),
    threads: 2 # TODO
    resources:
        mem="12GB",
        runtime="10m",
    shell:
        logged(
            "./{input.script}"
            "  {input.palmer_reads}"
            "  {input.ref_mei}"
            "  {input.orig_mei}"
            "  {wildcards.mei}"
            "  {output.out_dir}"
            "  {output.out_summary}"
        )


def pp_mei(wc):
    return {
        "LINE":  "meis/union/L1.inter.fi",
        "AluYa": "meis/union/ALU.inter.fi",
        "AluYb": "meis/union/ALU.inter.fi",
    }[wc.mei]

rule intersect_again:
    log:
        scratch("_logs/intersect_again/{sample}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/intersect_again/{sample}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        script="scripts/intersect-again.sh",
        container=containers("nanopal-binaries"),
        bam=scratch("alignment/{sample}/alignment.bam"),
        ref_mei=ref_mei,
        pp_mei=pp_mei,
        in_summary=scratch("intersect/{sample}/{mei}/summary.final.txt"),
        valid_read_ids=scratch("find_valid_reads/{sample}/RC.all.list"),
    output:
        out_dir=directory(scratch("intersect_again/{sample}/{mei}/")),
        out_summary=scratch("intersect_again/{sample}/{mei}/summary.final.2.txt"),
        out_result_log=scratch("intersect_again/{sample}/{mei}/result-log.txt"),
    threads: 2 # TODO
    resources:
        mem="12GB",
        runtime="10m",
    shell:
        logged(
            "./{input.script}"
            "  {input.bam}"
            "  {input.ref_mei}"
            "  {input.pp_mei}"
            "  {input.in_summary}"
            "  {input.valid_read_ids}"
            "  {wildcards.mei}"
            "  {output.out_dir}"
            "  {output.out_summary}"
            "  {output.out_result_log}"
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
        expand(scratch("alignment/{sample}/alignment.bam.bai"), sample=config["samples"]),

rule _palmer:
    localrule: True
    input:
        expand(
            scratch("palmer/{sample}/{mei}/{chromosome}/done"),
            sample=config["samples"],
            mei=config["mobile_elements"],
            chromosome=config["chromosomes"],
        ),

rule _on_target:
    localrule: True
    input:
        expand(
            scratch("find_on_target/{sample}/{mei}/read.all.txt"),
            sample=config["samples"],
            mei=config["mobile_elements"],
        ),
        expand(
            scratch("palmer_on_target/{sample}/{mei}/read.all.palmer.final.txt"),
            sample=config["samples"],
            mei=config["mobile_elements"],
        ),

rule _intersect:
    localrule: True
    input:
        expand(
            scratch("intersect/{sample}/{mei}/summary.final.txt"),
            sample=config["samples"],
            mei=config["mobile_elements"],
        ),
        expand(
            scratch("intersect_again/{sample}/{mei}/summary.final.2.txt"),
            sample=config["samples"],
            mei=config["mobile_elements"],
        ),

rule _all:
    localrule: True
    input:
        rules._palmer.input,
        rules._alignment.input,
        rules._on_target.input,
        rules._intersect.input,
