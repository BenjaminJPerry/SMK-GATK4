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


SAMPLES, = glob_wildcards("fastq/MFF/{samples}_R1.fastq.gz")


rule all:
    input:
        expand("results/01_mapping/{samples}.sorted.mkdups.bam", samples = SAMPLES),


rule bwa_mem:
    input:
        read1 = "fastq/MFF/{samples}_R1.fastq.gz",
        read2 = "fastq/MFF/{samples}_R2.fastq.gz",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna",
    output: 
        temp("results/01_mapping/{samples}.bam")
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
        partition="compute",
    shell:
        """
        ANIMAL="$(echo {wildcards.samples} | cut -d "_" -f 1)"
        #RUN="$(echo {wildcards.samples} | cut -d "_" -f 4)"
        #LANE="$(echo {wildcards.samples} | cut -d "_" -f 3)"
        READGROUP="@RG\\tID:$ANIMAL\\tPU:MFF\\tSM:$ANIMAL\\tPL:ILLUMINA\\tLB:ILLUMINA"

        bwa mem -Y -R $READGROUP -t {threads} -K 10000000 {input.referenceGenome} {input.read1} {input.read2} | samtools view --threads {threads} -bS -o {output}

        """


rule samtools_sort:
    input:
        bam = "results/01_mapping/{samples}.bam",
    output:
        sortedbam = temp("results/01_mapping/{samples}.sorted.bam"),
    log:
        "logs/samtools_sort_bwa.{samples}.log"
    benchmark:
        "benchmarks/samtools_sort_bwa.{samples}.tsv"
    conda:
        "bwa"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt -1) * 24),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute"
    shell:
        """

        samtools sort -l 8 -m 2G --threads {threads} {input.bam} > {output.sortedbam} && samtools index {output.sortedbam}
        
        """


rule gatk_MarkDuplicates:
    input:
        "results/01_mapping/{samples}.sorted.bam"
    output:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        metrics = "results/01_mapping/{samples}_mkdups_metrics.txt"
    log:
        "logs/gatk_MarkDuplicates.{samples}.log"
    benchmark:
        "benchmarks/gatk_MarkDuplicates.{samples}.tsv"
    conda:
        "gatk-4.5.0.0"
    threads:2
    resources:
        mem_gb = lambda wildcards, attempt: 128 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'gatk --java-options "-Xms2G -Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'MarkDuplicates '
        '-I {input} '
        '-O {output.bam} '
        '-M {output.metrics} '
        '--TMP_DIR {resources.DTMP} '
        '--COMPRESSION_LEVEL 8 '
        '&> {log}.attempt.{resources.attempt} '