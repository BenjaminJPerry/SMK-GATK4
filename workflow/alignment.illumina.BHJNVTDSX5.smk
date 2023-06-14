# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

#configfile: "config/config.yaml"


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


SAMPLES, = glob_wildcards("fastq/BHJNVTDSX5/{samples}_R1.fastq.gz")


rule all:
    input:
        expand("results/01_mapping/{samples}.sorted.bam", samples = SAMPLES),


rule bwa_mem:
    input:
        read1 = "fastq/BHJNVTDSX5/{samples}_R1.fastq.gz",
        read2 = "fastq/BHJNVTDSX5/{samples}_R2.fastq.gz",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output: 
        temp("results/01_mapping/{samples}.bam"),
    log:
        "logs/bwa_mem.{samples}.log"
    benchmark:
        "benchmarks/bwa_mem.{samples}.tsv"
    conda:
        "bwa"
    threads:24
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 24),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition="large,milan",
    shell:
        """
        ANIMAL="$(echo {wildcards.samples} | cut -d "_" -f 1)"
        RUN="$(echo {wildcards.samples} | cut -d "_" -f 4)"
        LANE="$(echo {wildcards.samples} | cut -d "_" -f 3)"
        READGROUP="@RG\\tID:$ANIMAL\\tPU:$RUN.$LANE\\tSM:$ANIMAL\\tPL:ILLUMINA\\tLB:ILLUMINA"

        bwa mem -Y -R $READGROUP -t {threads} -K 10000000 {input.referenceGenome} {input.read1} {input.read2} | samtools view --threads {threads} -bS -o {output}

        """


rule samtools_sort:
    input:
        bam = "results/01_mapping/{samples}.bam",
    output:
        sortedbam = "results/01_mapping/{samples}.sorted.bam",
    log:
        "logs/samtools_sort_bwa.{samples}.log"
    benchmark:
        "benchmarks/samtools_sort_bwa.{samples}.tsv"
    conda:
        "bwa"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempts: 26 + ((attempt -1) * 24),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan"
    shell:
        """

        samtools sort -l 8 -m 2G --threads {threads} {input.bam} > {output.bam} && samtools index {output.bam}

        """

