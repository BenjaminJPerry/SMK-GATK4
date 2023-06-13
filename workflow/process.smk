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


rule merge_bams:


rule samtools_sort: #TODO Make snakemake pipe
    input:
        "",
    output:
        sorted_bam = "results/01_mapping/{samples}.sorted.bam",
    log:
        "logs/samtools_sort.{samples}.log"
    conda:
        "samtools"
    threads: 8
    resources:
    shell:
        "gatk SortSam "
        "--java-options '-Xmx{resources.mem_gb}G' "
        "< {input.bam} "
        "-O {output.sorted_bam} "
        "2> {log}"


rule gatk_MarkDuplicates:
    input:
        "results/mapped/{sample}_sorted.bam"
    output:
        bam = temp("results/mapped/{sample}_sorted_mkdups.bam"),
        metrics = "results/mapped/{sample}_sorted_mkdups_metrics.txt"
    params:
        tdir = config["TEMPDIR"]
    log:
        "logs/gatk_MarkDuplicates.{sample}.log"
    benchmark:
        "benchmarks/gatk_MarkDuplicates/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 2
    resources:
        mem_gb = 8,
        # partition = config["PARTITION"]["CPU"]
    message:
        "Locating and tagging duplicate reads in {input}"
    shell:
        "gatk MarkDuplicates "
        "--java-options '-Xmx{resources.mem_gb}G' "
        "-I {input} "
        "-O {output.bam} "
        "-M {output.metrics} "
        "--TMP_DIR {params.tdir} "
        "&> {log}"


