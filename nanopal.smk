import os


# Config ---------------------------------------------------------------------


configfile: "config/base.json"


IDS = config["datasets"].keys()
BATCH_ID = config["batch_id"]
EXCLUSIONS = config["exclusions"]
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
def dataset_input_files(wc):
    paths = config["datasets"][wc.id]
    if isinstance(paths, str):
        paths = [paths]
    return paths

rule input:
    log:
        scratch("_logs/input/{id}.log"),
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        script="scripts/retrieve-input.sh",
        input_files=dataset_input_files,
    output:
        fastq=scratch("{id}/input/batch.fastq"),
        fasta=scratch("{id}/input/batch.fasta"),
    params:
        input_file_string=lambda wc, input: '\n'.join(input.input_files),
    threads: 2
    resources:
        mem="4GB",
        runtime="1h",
    shell:
        logged(
            "./{input.script} "
            "  {params.input_file_string:q}"
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

rule detect_ligation_artifacts:
    log:
        scratch("_logs/detect_ligation_artifacts/{id}.log"),
    benchmark:
        scratch("_benchmarks/detect_ligation_artifacts/{id}.tsv")
    container:
        containers("liger2liger")
    input:
        container=containers("liger2liger"),
        fastq=scratch("{id}/input/batch.fastq"),
        index=scratch("reference/index.mmi"),
    output:
        result=scratch("{id}/detect_ligation_artifacts/ligation_artifacts.txt"),
    params:
        output_dir=scratch("{id}/detect_ligation_artifacts"),
    threads: 24
    resources:
        mem="72GB",
        runtime="3h",
    shell:
        logged(
            "cd {params.output_dir}",
            "liger2liger --ref {input.index} --fastq {input.fastq} --n_threads {threads}",
            "cp {params.output_dir}/batch_VS_index/batch_VS_index.chimeric_reads.txt {output.result}"
        )

rule detect_hallucination:
    log:
        scratch("_logs/detect_hallucination/{id}.log"),
    benchmark:
        scratch("_benchmarks/detect_hallucination/{id}.tsv")
    container:
        containers("delulu")
    input:
        container=containers("delulu"),
        fastq=scratch("{id}/input/batch.fastq"),
    output:
        result=scratch("{id}/detect_hallucination/hallucination.csv"),
    threads: 6
    resources:
        mem="12GB",
        runtime="1h",
    shell:
        logged(
            "delulu --threads {threads} {input.fastq} > {output.result}"
        )

rule dump_alignments_for_minimera:
    log:
        scratch("_logs/dump_alignments_for_minimera/{id}.log"),
    benchmark:
        scratch("_benchmarks/dump_alignments_for_minimera/{id}.tsv")
    container:
        containers("minimap2") # Hack because it has samtools and gawk
    input:
        container=containers("minimap2"),
        script="scripts/dump-alignments-for-minimera.sh",
        bam=scratch("{id}/alignment/alignment.bam"),
    output:
        result=scratch("{id}/dump_alignments_for_minimera/alignments.csv"),
    threads: 1
    resources:
        mem="4GB",
        runtime="1h",
    shell:
        logged(
            "./{input.script} "
            "  {input.bam}"
            "> {output.result}"
        )

rule minimera:
    log:
        scratch("_logs/minimera/{id}.log"),
    benchmark:
        scratch("_benchmarks/minimera/{id}.tsv")
    container:
        containers("minimera")
    input:
        container=containers("minimera"),
        fastq=scratch("{id}/input/batch.fastq"),
        # alignments=scratch("{id}/dump_alignments_for_minimera/alignments.csv"),
    output:
        result_csv=scratch("{id}/minimera/foldbacks.csv"),
        # result_bed=scratch("{id}/minimera/foldbacks.bed"),
    params:
        output_dir=scratch("{id}/minimera"),
    threads: 18
    resources:
        mem="48GB",
        runtime="4h",
    shell:
        logged(
            "minimera "
            "  {input.fastq}"
            "  --threads {threads}"
            "  --output {params.output_dir}"
            "  --no-progress"
            "  --no-plot-foldbacks"
            "  --monotony-threshold 0.5"
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
        fastq=scratch("{id}/input/batch.fastq"),
        index=scratch("reference/index.mmi"),
    output:
        bam=scratch("{id}/alignment/alignment.bam"),
        bai=scratch("{id}/alignment/alignment.bam.bai"),
    params:
        output_dir=scratch("{id}/alignment"),
        samtools_threads=lambda wc, threads: min(16, threads),
    threads: 24
    resources:
        mem="64GB",
        runtime="4h",
    shell:
        logged(
            "rm -f {output.bam}",
            "find {params.output_dir} -name '*.tmp.*' | xargs -r rm",
            "minimap2"
            "  -a"
            "  -y"
            "  -x map-ont"
            "  -t {threads}"
            "  --secondary=no"
            "  --eqx"
            "  -Y"
            "  {input.index} {input.fastq}"
            "  | samtools sort -o {output.bam} -@ {params.samtools_threads}",
            "samtools index {output.bam} -@ {params.samtools_threads}"
        )

rule find_revcomp_reads:
    log:
        scratch("_logs/find_revcomp_reads/{id}.log"),
    benchmark:
        scratch("_benchmarks/find_revcomp_reads/{id}.tsv")
    container:
        containers("samtools")
    input:
        container=containers("samtools"),
        bam=scratch("{id}/alignment/alignment.bam"),
    output:
        revcomp_read_ids=scratch("{id}/find_revcomp_reads/RC.all.list"),
    threads: 2
    resources:
        mem="2GB",
        runtime="30m",
    shell:
        logged(
            "samtools view {input.bam}"
            "    -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 -f 0x10"
            "  | awk '{{print $1}}'"
            "  > {output.revcomp_read_ids}"
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
        bam=scratch("{id}/alignment/alignment.bam"),
    output:
        blast_results=scratch("{id}/{mei}/palmer/{chromosome}/collected_blastn_refine.txt")
    params:
        workdir=scratch("{id}/{mei}/palmer/{chromosome}/"),
        reference_version=config["reference_version"],
        palmer_mei=palmer_mei_param,
    threads: 2  # TODO
    resources:
        mem="8GB",
        runtime="4h",
    shell:
        # First run Palmer.  Then grab all the Palmer blast results from the
        # regional subsets in the working directories.  Finally, clean up the
        # Palmer directory once we've gotten what we need from it to avoid
        # exhausting all the inodes on the drive when running many datasets.
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
            " && find {params.workdir} -mindepth 1 -name blastn_refine.txt"
            "    | xargs cat > {output.blast_results}"
            " && find {params.workdir} -mindepth 1 -maxdepth 2 -name 'chr*_*_*' -type d"
            "    | xargs rm -r"
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
        palmer_blasts=expand(
            scratch("{{id}}/{{mei}}/palmer/{chromosome}/collected_blastn_refine.txt"),
            chromosome=config["chromosomes"],
        ),
        bam=scratch("{id}/alignment/alignment.bam"),
    output:
        palmer_blast=scratch("{id}/{mei}/gather_matches/blastn_refine.all.txt"),
        palmer_cigar=scratch("{id}/{mei}/gather_matches/mapped.info.txt"),
    params:
        palmer_dir=scratch("{id}/{mei}/palmer/"),
    threads: 1
    resources:
        mem="2GB",
    shell:
        logged(
            "./{input.script} "
            "  {output.palmer_blast}"
            "  {output.palmer_cigar}"
            "  {input.bam}"
            "  {input.palmer_blasts}"
        )

rule parse_cigar:
    log:
        scratch("_logs/parse_cigar/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/parse_cigar/{id}_{mei}.tsv")
    container:
        containers("nanopal-binaries")
    input:
        container=containers("nanopal-binaries"),
        cigar_matches=scratch("{id}/{mei}/gather_matches/mapped.info.txt"),
    output:
        cigar_results=scratch("{id}/{mei}/parse_cigar/cigar_results.all.txt"),
        cigar_ref=scratch("{id}/{mei}/parse_cigar/cigar_ref.txt"),
        mapped_info=scratch("{id}/{mei}/parse_cigar/mapped.info.final.txt"),
    threads: 2
    resources:
        mem="2GB",
        runtime="30m",
    shell:
        logged(
            "awk '{{print $4}}' {input.cigar_matches} | cigar_parser > {output.cigar_results}",
            "awk '{{print $4+$6+$10}}' {output.cigar_results} > {output.cigar_ref}",
            "paste {input.cigar_matches} {output.cigar_ref}"
            "  | awk '{{print $1, $2, $3, $3+$5}}'"
            "  > {output.mapped_info}"
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
        container=containers("blast"),
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
        reads_fasta=scratch("{id}/input/batch.fasta"),
    output:
        nanopal_reads=scratch("{id}/{mei}/find_on_target/read.all.txt")
    params:
        mei_cut_site=mei_cut_site,
        blast_threads=lambda wc, threads: max(threads-2, 1),
    threads: 12
    resources:
        mem="24GB",
        runtime="4h",
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
    input:
        script="scripts/palmer-on-target.sh",
        palmer_blast=scratch("{id}/{mei}/gather_matches/blastn_refine.all.txt"),
        palmer_map=scratch("{id}/{mei}/parse_cigar/mapped.info.final.txt"),
        nanopal_reads=scratch("{id}/{mei}/find_on_target/read.all.txt"),
    output:
        scratch("{id}/{mei}/palmer_on_target/read.all.palmer.final.txt")
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
    if EXCLUSIONS == "none":
        return "meis/empty.txt"
    else:
        return {
            "GM12878": {
                "LINE":  "meis/PALMER.NA12878.L1.txt",
                "AluYa": "meis/PALMER.NA12878.ALU.txt",
                "AluYb": "meis/PALMER.NA12878.ALU.txt",
                "SVA_E": "meis/PALMER.NA12878.SVA.txt",
                "SVA_F": "meis/PALMER.NA12878.SVA.txt",
            }[wc.mei]
        }[EXCLUSIONS]

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
        palmer_reads=scratch("{id}/{mei}/palmer_on_target/read.all.palmer.final.txt"),
        ref_mei=ref_mei,
        orig_mei=orig_mei,
    output:
        out_dir=directory(scratch("{id}/{mei}/intersect/")),
        out_summary=scratch("{id}/{mei}/intersect/summary.final.txt"),
    threads: 2 # TODO
    resources:
        mem="27GB",
        runtime="30m",
    shell:
        logged(
            "./{input.script}"
            "  {input.palmer_reads}"
            "  {input.ref_mei}"
            "  {input.orig_mei}"
            "  {output.out_dir}"
            "  {output.out_summary}"
        )

def pp_mei(wc):
    if EXCLUSIONS == "none":
        return "meis/union/empty.txt"
    else:
        return {
            "GM12878": {
                "LINE":  "meis/union/L1.inter.fi",
                "AluYa": "meis/union/ALU.inter.fi",
                "AluYb": "meis/union/ALU.inter.fi",
                "SVA_E": "meis/union/SVA.inter.fi",
                "SVA_F": "meis/union/SVA.inter.fi",
            }[wc.mei]
        }[EXCLUSIONS]

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
        in_summary=scratch("{id}/{mei}/intersect/summary.final.txt"),
        revcomp_read_ids=scratch("{id}/find_revcomp_reads/RC.all.list"),
        foldbacks=scratch("{id}/minimera/foldbacks.csv"),
        ligation_artifacts=scratch("{id}/detect_ligation_artifacts/ligation_artifacts.txt"),
        hallucination=scratch("{id}/detect_hallucination/hallucination.csv"),
    output:
        out_dir=directory(scratch("{id}/{mei}/intersect_again/")),
        out_summary=scratch("{id}/{mei}/intersect_again/summary.final.2.txt"),
        out_result_log=scratch("{id}/{mei}/intersect_again/result-log.txt"),
        out_potential_meis=scratch("{id}/{mei}/intersect_again/potential.clustered.txt.fi"),
    threads: 2 # TODO
    resources:
        mem="12GB",
        runtime="30m",
    shell:
        logged(
            "./{input.script}"
            "  {wildcards.id}"
            "  {input.ref_mei}"
            "  {input.pp_mei}"
            "  {input.in_summary}"
            "  {input.revcomp_read_ids}"
            "  {input.foldbacks}"
            "  {input.ligation_artifacts}"
            "  {input.hallucination}"
            "  {wildcards.mei}"
            "  {output.out_dir}"
            "  {output.out_summary}"
            "  {output.out_result_log}"
        )

rule output_results:
    log:
        scratch("_logs/output_results/{id}_{mei}.log"),
    benchmark:
        scratch("_benchmarks/output_results/{id}_{mei}.tsv")
    container:
        containers("samtools")
    input:
        script="scripts/output-results.sh",
        container=containers("nanopal-binaries"),
        result_log=scratch("{id}/{mei}/intersect_again/result-log.txt"),
        potential_meis=scratch("{id}/{mei}/intersect_again/potential.clustered.txt.fi"),
        bam=scratch("{id}/alignment/alignment.bam"),
    output:
        csv=scratch("{id}/{mei}/output_results/result.csv"),
        log=scratch("{id}/{mei}/output_results/result.txt"),
        bed=scratch("{id}/{mei}/output_results/nanopal_calls-all-{id}-{mei}.bed"),
        bed_multi=scratch("{id}/{mei}/output_results/nanopal_calls-multiple_support-{id}-{mei}.bed"),
    params:
        out_dir=scratch("{id}/{mei}/output_results/")
    threads: 16
    resources:
        mem="32GB",
        runtime="3h",
    shell:
        logged(
            "./{input.script}"
            "  {wildcards.id}"
            "  {wildcards.mei}"
            "  {input.bam}"
            "  {input.potential_meis}"
            "  {input.result_log}"
            "  {threads}"
            "  {params.out_dir}"
            "  {output.csv}"
            "  {output.log}"
            "  {output.bed}"
            "  {output.bed_multi}"
        )

rule collect_results:
    localrule: True
    log:
        scratch("_logs/collect_results.log"),
    benchmark:
        scratch("_benchmarks/collect_results.tsv")
    input:
        script="scripts/collect-results.sh",
        result_csvs=expand(
            scratch("{id}/{mei}/output_results/result.csv"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
        result_beds=expand(
            scratch("{id}/{mei}/output_results/nanopal_calls-all-{id}-{mei}.bed"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
    output:
        out_csv=scratch("collect-results/results.csv"),
    params:
        out_dir=scratch("collect-results/"),
        result_dirs=expand(
            scratch("{id}/{mei}/output_results/"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
    threads: 1
    shell:
        logged(
            "./{input.script}"
            "  {params.out_dir}"
            "  {output.out_csv}"
            "  {params.result_dirs}"
        )


# PHONY -----------------------------------------------------------------------
rule _input:
    localrule: True
    input:
        expand(scratch("{id}/input/batch.fastq"), id=IDS),
        expand(scratch("{id}/input/batch.fasta"), id=IDS),

rule _alignment:
    localrule: True
    input:
        expand(scratch("{id}/alignment/alignment.bam"), id=IDS),
        expand(scratch("{id}/alignment/alignment.bam.bai"), id=IDS),

rule _palmer:
    localrule: True
    input:
        expand(
            scratch("{id}/{mei}/palmer/{chromosome}/collected_blastn_refine.txt"),
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
            scratch("{id}/{mei}/find_on_target/read.all.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
        expand(
            scratch("{id}/{mei}/palmer_on_target/read.all.palmer.final.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),

rule _intersect:
    localrule: True
    input:
        expand(
            scratch("{id}/{mei}/intersect/summary.final.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),
        expand(
            scratch("{id}/{mei}/intersect_again/summary.final.2.txt"),
            id=IDS,
            mei=config["mobile_elements"],
        ),

rule _ligation_artifacts:
    localrule: True
    input:
        expand(
            scratch("{id}/detect_ligation_artifacts/ligation_artifacts.txt"),
            id=IDS,
        ),

rule _hallucination:
    localrule: True
    input:
        expand(
            scratch("{id}/detect_hallucination/hallucination.txt"),
            id=IDS,
        ),

rule _minimera:
    localrule: True
    input:
        expand(scratch("{id}/minimera/foldbacks.csv"), id=IDS),
        # expand(scratch("{id}/minimera/foldbacks.bed"), id=IDS),

rule _results:
    localrule: True
    input:
        scratch("collect-results/results.csv"),
        rules._minimera.input,
        rules._ligation_artifacts.input,

rule _all:
    default_target: True
    localrule: True
    input:
        rules._results.input,
