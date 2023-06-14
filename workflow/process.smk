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


SAMPLES, = glob_wildcards("results/01_mapping/{samples}.bam")


rule all:
    input:
        "",


rule samtools_merge:
    input:
        bam = "fastq/MGI/{samples}_R1.fastq.gz",
    output: 
        "results/01_mapping/{samples}.merged.bam"
    log:
        "logs/samtools_merge.{samples}.log"
    benchmark:
        "benchmarks/samtools_merge.{samples}.tsv"
    conda:
        "bwa"
    threads:24
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 24),
        time = lambda wildcards, attempt: 360 + ((attempt - 1) * 60),
        partition="large,milan",
    shell:
        """
        ANIMAL="$(echo {wildcards.samples} | cut -d "_" -f 1)"
        RUN="$(echo {wildcards.samples} | cut -d "_" -f 2)"
        LANE="$(echo {wildcards.samples} | cut -d "_" -f 3)"
        READGROUP="@RG\\tID:$ANIMAL\\tPU:$RUN.$LANE\\tSM:$ANIMAL\\tPL:DNBSEQ\\tLB:DNBSEQ"

        bwa mem -Y -R $READGROUP -t {threads} -K 100000000 {input.referenceGenome} {input.read1} {input.read2} | samtools view --threads {threads} -bS -o {output}

        """


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


