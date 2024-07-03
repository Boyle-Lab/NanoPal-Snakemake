import os


# Config ---------------------------------------------------------------------


configfile: "config/base.json"


IDS = config["datasets"].keys()
BATCH_ID = config["batch_id"]
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
    commands = "\n".join(('date',) + shell_commands + ('date',))
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

def dataset_input_dirs(wc):
    return [directory(d) for d in config["datasets"][wc.id]]

rule input:
    log:
        scratch("_logs/input/{id}.log"),
    input:
        script="scripts/retrieve-input.sh",
        input_dirs=dataset_input_dirs,
    output:
        fastq=scratch("input/{id}/batch.fastq"),
        fasta=scratch("input/{id}/batch.fasta"),
        # basecall_info=scratch("input/{id}/basecall_info.txt"),
    params:
        input_dir_string=lambda wc, input: '\n'.join(input.input_dirs),
    threads: 2
    resources:
        mem="4GB",
        runtime="1h",
    shell:
        logged(
            "./{input.script} "
            "  {params.input_dir_string:q}"
            "  {output.fastq}"
            "  {output.fasta}"
            # "  {output.basecall_info}"
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
        scratch("_logs/alignment/{id}.log"),
    benchmark:
        scratch("_benchmarks/alignment/{id}.tsv")
    container:
        containers("minimap2")
    input:
        container=containers("minimap2"),
        fastq=scratch("input/{id}/batch.fastq"),
        index=scratch("reference/index.mmi"),
    output:
        bam=scratch("alignment/{id}/alignment.bam"),
        bai=scratch("alignment/{id}/alignment.bam.bai"),
    params:
        output_dir=scratch("alignment/{id}"),
        samtools_threads=lambda wc, threads: min(16, threads),
    threads: 36
    resources:
        mem="128GB",
        runtime="3h",
    shell:
        logged(
            "rm -f {output.bam}",
            "find {params.output_dir} -name '*.tmp.*' | xargs -r rm",
            "minimap2"
            "  -a"
            "  -x map-ont"
            "  -t {threads}"
            "  --secondary=no"
            "  --eqx"
            "  -Y"
            "  {input.index} {input.fastq}"
            "  | samtools sort -o {output.bam} -@ {params.samtools_threads}",
            "samtools index {output.bam} -@ {params.samtools_threads}"
        )

rule find_valid_reads:
    localrule: True
    log:
        scratch("_logs/find_valid_reads/{id}.log"),
    benchmark:
        scratch("_benchmarks/find_valid_reads/{id}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        bam=scratch("alignment/{id}/alignment.bam"),
    output:
        valid_read_ids=scratch("find_valid_reads/{id}/RC.all.list"),
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
        'SVA_E': 'SVA',
        'SVA_F': 'SVA',
    }[wc.mei]

rule palmer:
    log:
        scratch("_logs/palmer/{id}_{mei}_{chromosome}.log"),
    benchmark:
        scratch("_benchmarks/palmer/{id}_{mei}_{chromosome}.tsv")
    container:
        containers("palmer")
    input:
        container=containers("palmer"),
        reference=config["reference"],
        bam=scratch("alignment/{id}/alignment.bam"),
    output:
        sentinel=touch(scratch("palmer/{id}/{mei}/{chromosome}/done")),
    params:
        workdir=scratch("palmer/{id}/{mei}/{chromosome}/"),
        reference_version=config["reference_version"],
        palmer_mei=palmer_mei_param,
    threads: 2  # TODO
    resources:
        mem="8GB",
        runtime="4h",
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
        scratch("_logs/gather_matches/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/gather_matches/{id}_{mei}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        script="scripts/gather-palmer-results.sh",
        palmer_sentinels=expand(
            scratch("palmer/{{id}}/{{mei}}/{chromosome}/done"),
            chromosome=config["chromosomes"],
        ),
        bam=scratch("alignment/{id}/alignment.bam"),
    output:
        palmer_blast=scratch("gather_matches/{id}/{mei}/blastn_refine.all.txt"),
        palmer_cigar=scratch("gather_matches/{id}/{mei}/mapped.info.txt"),
    params:
        palmer_dir=scratch("palmer/{id}/{mei}/"),
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
        scratch("_logs/parse_cigar/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/parse_cigar/{id}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        container=containers("nanopal-binaries"),
        cigar_matches=scratch("gather_matches/{id}/{mei}/mapped.info.txt"),
    output:
        cigar_results=scratch("parse_cigar/{id}/{mei}/cigar_results.all.txt"),
        cigar_ref=scratch("parse_cigar/{id}/{mei}/cigar_ref.txt"),
        mapped_info=scratch("parse_cigar/{id}/{mei}/mapped.info.final.txt"),
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
        'SVA_E': "meis/SVA_E",
        'SVA_F': "meis/SVA_F",
    }[wc.mei]

rule mei_db:
    localrule: True
    log:
        scratch("_logs/mei_db/{mei}.log")
    container:
        containers("blast")
    input:
        container=containers("nanopal-binaries"),
        mei_fasta=mei_fasta,
    output:
        mei_db=multiext(scratch("mei_db/{mei}/{mei}"), ".ndb", ".nhr", ".nin", ".njs", ".not", ".nsq", ".ntf", ".nto"),
        mei_marker=touch(scratch("mei_db/{mei}/{mei}")),
    params:
        output_prefix=scratch("mei_db/{mei}/{mei}")
    threads: 1
    shell:
        logged(
            "makeblastdb"
            "  -in {input.mei_fasta}"
            "  -dbtype nucl"
            "  -out {params.output_prefix}"
        )

def mei_cut_site(wc):
    return {
        'LINE':  5900,
        'AluYa':  225,
        'AluYb':  227,
        'SVA_E':  930,
        'SVA_F': 1006,
    }[wc.mei]

rule find_on_target:
    log:
        scratch("_logs/find_on_target/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/find_on_target/{id}_{mei}.tsv")
    container:
        containers("blast")
    input:
        container=containers("blast"),
        script="scripts/blast-reads.sh",
        mei_db=scratch("mei_db/{mei}/{mei}"),
        reads_fasta=scratch("input/{id}/batch.fasta"),
    output:
        nanopal_reads=scratch("find_on_target/{id}/{mei}/read.all.txt")
    params:
        mei_cut_site=mei_cut_site,
        blast_threads=lambda wc, threads: max(threads-2, 1),
    threads: 10
    resources:
        mem="24GB",
        runtime="1h",
    shell:
        logged(
            "./{input.script}"
            "  {input.mei_db}"
            "  {input.reads_fasta}"
            "  {params.mei_cut_site}"
            "  {params.blast_threads}"
            "  > {output}"
        )

rule palmer_on_target:
    localrule: True
    log:
        scratch("_logs/palmer_on_target/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/palmer_on_target/{id}_{mei}.tsv")
    container:
        containers("palmer")
    input:
        container=containers("palmer"),
        script="scripts/palmer-on-target.sh",
        palmer_blast=scratch("gather_matches/{id}/{mei}/blastn_refine.all.txt"),
        palmer_map=scratch("parse_cigar/{id}/{mei}/mapped.info.final.txt"),
        nanopal_reads=scratch("find_on_target/{id}/{mei}/read.all.txt"),
    output:
        scratch("palmer_on_target/{id}/{mei}/read.all.palmer.final.txt")
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
        "SVA_E": "meis/hg38.RM.SVA.ref",
        "SVA_F": "meis/hg38.RM.SVA.ref",
    }[wc.mei]

def orig_mei(wc):
    return {
        "LINE":  "meis/PALMER.NA12878.L1.txt",
        "AluYa": "meis/PALMER.NA12878.ALU.txt",
        "AluYb": "meis/PALMER.NA12878.ALU.txt",
        "SVA_E": "meis/PALMER.NA12878.SVA.txt",
        "SVA_F": "meis/PALMER.NA12878.SVA.txt",
    }[wc.mei]

rule intersect:
    log:
        scratch("_logs/intersect/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/intersect/{id}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        script="scripts/intersect.sh",
        container=containers("nanopal-binaries"),
        palmer_reads=scratch("palmer_on_target/{id}/{mei}/read.all.palmer.final.txt"),
        ref_mei=ref_mei,
        orig_mei=orig_mei,
    output:
        out_dir=directory(scratch("intersect/{id}/{mei}/")),
        out_summary=scratch("intersect/{id}/{mei}/summary.final.txt"),
    threads: 2 # TODO
    resources:
        mem="12GB",
        runtime="30m",
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
        "SVA_E": "meis/union/SVA.inter.fi",
        "SVA_F": "meis/union/SVA.inter.fi",
    }[wc.mei]

rule intersect_again:
    log:
        scratch("_logs/intersect_again/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/intersect_again/{id}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        script="scripts/intersect-again.sh",
        container=containers("nanopal-binaries"),
        ref_mei=ref_mei,
        pp_mei=pp_mei,
        in_summary=scratch("intersect/{id}/{mei}/summary.final.txt"),
        valid_read_ids=scratch("find_valid_reads/{id}/RC.all.list"),
    output:
        out_dir=directory(scratch("intersect_again/{id}/{mei}/")),
        out_summary=scratch("intersect_again/{id}/{mei}/summary.final.2.txt"),
        out_result_log=scratch("intersect_again/{id}/{mei}/result-log.txt"),
    threads: 2 # TODO
    resources:
        mem="12GB",
        runtime="30m",
    shell:
        logged(
            "./{input.script}"
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
        expand(scratch("input/{id}/batch.fastq"), id=IDS),
        expand(scratch("input/{id}/batch.fasta"), id=IDS),

rule _alignment:
    localrule: True
    input:
        expand(scratch("alignment/{id}/alignment.bam"), id=IDS),
        expand(scratch("alignment/{id}/alignment.bam.bai"), id=IDS),

rule _palmer:
    localrule: True
    input:
        expand(
            scratch("palmer/{id}/{mei}/{chromosome}/done"),
            id=IDS,
            mei=config["mobile_elements"],
            chromosome=config["chromosomes"],
        ),

rule _mei_db:
    localrule: True
    input:
        expand(
            scratch("mei_db/{mei}/{mei}"),
            mei=config["mobile_elements"]
        ),

rule _on_target:
    localrule: True
    input:
        expand(
            scratch("find_on_target/{id}/{mei}/read.all.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
        expand(
            scratch("palmer_on_target/{id}/{mei}/read.all.palmer.final.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),

rule _intersect:
    localrule: True
    input:
        expand(
            scratch("intersect/{id}/{mei}/summary.final.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
        expand(
            scratch("intersect_again/{id}/{mei}/summary.final.2.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),

rule _all:
    default_target: True
    localrule: True
    input:
        rules._palmer.input,
        rules._alignment.input,
        rules._on_target.input,
        rules._intersect.input,
