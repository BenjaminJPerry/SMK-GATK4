# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"


import os


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


SAMPLES, = glob_wildcards("fastq/AHJNKHDSX5{samples}_R1.fastq.gz")


rule all:
    input:
        "",

rule gatk_HaplotypeCaller_first_pass:
    input:
        bams = "results/mapped/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config["REFGENOME"]),
    output:
        protected("results/called/{sample}_raw_snps_indels.vcf")
    params:
        tdir = config["TEMPDIR"],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command
    log:
        "logs/gatk_HaplotypeCaller_single/{sample}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_single/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 2
    resources:
        mem_gb = 8,
        # partition = config["PARTITION"]["CPU"]
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bams}"
    shell:
        "gatk HaplotypeCaller "
        "--java-options '-Xmx{resources.mem_gb}G' "
        "-I {input.bams} "
        "-R {input.refgenome} "
        # "-D {input.dbsnp} "
        "-O {output} "
        "--tmp-dir {params.tdir} {params.padding} {params.intervals} "
        "&> {log}"


rule gatk_BaseRecalibrator:
    input:
        bams = "results/mapped/{sample}_sorted_mkdups.bam",
        refgenome = expand("{refgenome}", refgenome = config["REFGENOME"])
    output:
        report("results/mapped/{sample}_recalibration_report.grp", caption = "report/recalibration.rst", category = "Base recalibration")
    params:
        tdir = config["TEMPDIR"],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command,
        recalibration_resources = get_recal_resources_command
    log:
        "logs/gatk_BaseRecalibrator.{sample}.log"
    benchmark:
        "benchmarks/gatk_BaseRecalibrator/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 2
    resources:
        mem_gb = 8,
        # partition = config["PARTITION"]["CPU"]
    message:
        "Generating a recalibration table for {input.bams}"
    shell:
        "gatk BaseRecalibrator "
        "--java-options '-Xmx{resources.mem_gb}G' "
        "-I {input.bams} "
        "-R {input.refgenome} "
        "-O {output} "
        "--tmp-dir {params.tdir} {params.padding} {params.intervals} {params.recalibration_resources} "
        "&> {log}"


rule gatk_ApplyBQSR:
    input:
        bam = "results/mapped/{sample}_sorted_mkdups.bam",
        recal = "results/mapped/{sample}_recalibration_report.grp",
        refgenome = expand("{refgenome}", refgenome = config["REFGENOME"])
    output:
        bam = protected("results/mapped/{sample}_recalibrated.bam")
    params:
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command
    log:
        "logs/gatk_ApplyBQSR.{sample}.log"
    benchmark:
        "benchmarks/gatk_ApplyBQSR/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 2
    resources:
        mem_gb = 8,
        # partition = config["PARTITION"]["CPU"]
    message:
        "Applying base quality score recalibration and producing a recalibrated BAM file for {input.bam}"
    shell:
        "gatk ApplyBQSR "
        "--java-options '-Xmx{resources.mem_gb}G' "
        "-I {input.bam} "
        "-bqsr {input.recal} "
        "-R {input.refgenome} "
        "-O {output} {params.padding} {params.intervals}"


