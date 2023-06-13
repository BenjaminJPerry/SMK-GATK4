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


rule bwa_mem:
    input:
        read1 = "",
        read2 = "",
        referenceGenome = "",
    output: 
        "results/01_mapping/{samples}.bam"
    params:
        readgroup = 
    log:
        "logs/bwa_mem.{sample}.log"
    benchmark:
        "benchmarks/bwa_mem.{sample}.tsv"
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
        RUN="$(echo {wildcards.samples} | cut -d "_" -f 4)"
        LANE="$(echo {wildcards.samples} | cut -d "_" -f 3)"
        READGROUP="@RG\tID:$ANIMAL\tPU:$RUN.$LANE\tSM:$ANIMAL\tPL:ILLUMINA\tLB:ILLUMINA"

        bwa mem -Y -R $READGROUP -t {threads} -K 10000000 {input.referenceGenome} {input.read1} {input.read2} \
        | samtools view --threads {threads} -bS -o {output}

        """


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


