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

#TODO: Write a clever input_function to return the basenames for samples; platform level information is in the readgroup information
#Removed NEGATIVE,Blank due to low read mapping [GOOD]
SAMPLES = "120001,2201,2202,2203,2204,2205,2206,2207,2209,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,2223,2224,2225,2226,2227,849,864,Debonair,PDGY-16-16,PDGY-17-104,PDGY-18-17"


rule all:
    input:
#        expand("results/01_mapping/{samples}.merged.bam", samples = SAMPLES.split(","))
        expand("results/02_snvs/{samples}.raw.snvs.gvcf.gz", samples = SAMPLES.split(",")),


rule samtools_merge:
    priority: 100
    output: 
        temp("results/01_mapping/{samples}.merged.bam")
    log:
        "logs/samtools_merge.{samples}.log"
    benchmark:
        "benchmarks/samtools_merge.{samples}.tsv"
    conda:
        "bwa"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 128 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 720 + ((attempt - 1) * 360),
        partition="large,milan",
        tmpdir="temp"
    shell:
        """

        samtools merge -f --threads {threads} - results/01_mapping/{wildcards.samples}*.sorted.bam 2> {log} | samtools sort -l 8 -m 2G --threads {threads} > {output} && for file in $(ls results/01_mapping/{wildcards.samples}*.sorted.bam); do rm $file; done;
        
        for file in $(ls results/01_mapping/{wildcards.samples}*.sorted.bam.bai); do rm $file; done;
        
        samtools index {output}

        """

rule gatk_MarkDuplicates:
    input:
        "results/01_mapping/{samples}.merged.bam"
    output:
        bam = "results/01_mapping/{samples}.merged.sorted.mkdups.bam",
        metrics = "results/01_mapping/{samples}_sorted_mkdups_metrics.txt"
    log:
        "logs/gatk_MarkDuplicates.{samples}.log"
    benchmark:
        "benchmarks/gatk_MarkDuplicates.{samples}.tsv"
    threads:2
    resources:
        mem_gb = lambda wildcards, attempt: 128 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition="large,milan",
        tmpdir="/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.3.0.0-gimkl-2022a; '
        'gatk --java-options "-Xms{resources.mem_gb}G -Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'MarkDuplicates '
        '-I {input} '
        '-O {output.bam} '
        '-M {output.metrics} '
        '--TMP_DIR {resources.tmpdir} '
        '2>&1 {log}.attempt.{resources.attempt} '


rule gatk_HaplotypeCaller:
    input:
        bam = "results/01_mapping/{samples}.merged.sorted.mkdups.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        gvcf = "results/02_snvs/{samples}.raw.snvs.gvcf.gz",
    log:
        "logs/gatk_HaplotypeCaller_cohort.{samples}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller.{samples}.tsv"
    threads:2
    resources:
        mem_gb = lambda wildcards, attempt: 128 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition="large,milan",
        tmpdir="/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.3.0.0-gimkl-2022a ; '
        'gatk --java-options "-Xms{resources.mem_gb}G -Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'HaplotypeCaller '
        '-I {input.bam} '
        '-R {input.referenceGenome} '
        '-O {output.gvcf} '
        '-ERC GVCF '
        '--tmp-dir {resources.tmpdir} '
        '2>&1 {log}.attempt.{resources.attempt} '

