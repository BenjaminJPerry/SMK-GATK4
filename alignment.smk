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


FIDs, = glob_wildcards("fastq/{sample}.fastq.gz")


rule all:
    input:
        "",
        "",
        "",


rule bwa_mem:
    input:
        fastq = get_input_fastq,
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output: 
        temp("results/mapped/{sample}_sorted.bam")
    params:
        readgroup = "'@RG\\tID:{sample}_rg1\\tLB:lib1\\tPL:bar\\tSM:{sample}\\tPU:{sample}_rg1'",
        sortsam = "--MAX_RECORDS_IN_RAM=5000000 --SORT_ORDER=coordinate -I=/dev/stdin",
        tdir = expand("{tdir}", tdir = config['TEMPDIR'])
    log:
        "logs/bwa_mem.{sample}.log"
    benchmark:
        "benchmarks/bwa_mem.{sample}.tsv"
    conda:
        "bwa"
    threads: 12
    resources:
        mem_gb = 16,
        time = "08:00:00"
        # partition = config['PARTITION']['CPU']
    message:
        "Mapping sequences against a reference human genome with BWA-MEM for {input.fastq}"
    shell:
        'bwa mem '
        '-R {params.readgroup} '
        '-t {threads} '
        '-K 10000000 {input.refgenome} {input.fastq} | '
        'gatk SortSam '
        '--java-options "-Xmx{resources.mem_gb}G" '
        '{params.sortsam} '
        '-O={output} '
        '--TMP_DIR={params.tdir} '
        '&> {log}'


rule sort_bam: #TODO Make snakemake pipe
    input:
        "",
        "",
    output:
        "",
        "",
    log: ""
    conda: ""
    threads:
    resources:
        "",
        "",
    shell:
        "gatk ",
        "--java-options '-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_sort}G' ",
        "SortSam ",
        "",


rule gatk_MarkDuplicates:
    input:
        "results/mapped/{sample}_sorted.bam"
    output:
        bam = temp("results/mapped/{sample}_sorted_mkdups.bam"),
        metrics = "results/mapped/{sample}_sorted_mkdups_metrics.txt"
    params:
        tdir = config['TEMPDIR']
    log:
        "logs/gatk_MarkDuplicates.{sample}.log"
    benchmark:
        "benchmarks/gatk_MarkDuplicates/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 2
    resources:
        mem_gb = 8,
        # partition = config['PARTITION']['CPU']
    message:
        "Locating and tagging duplicate reads in {input}"
    shell:
        'gatk MarkDuplicates '
        '--java-options "-Xmx{resources.mem_gb}G" '
        '-I {input} '
        '-O {output.bam} '
        '-M {output.metrics} '
        '--TMP_DIR {params.tdir} '
        '&> {log}'


rule gatk_HaplotypeCaller_first_pass:
    input:
        bams = "results/mapped/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME']),
        # dbsnp = expand("{dbsnp}", dbsnp = config['dbSNP'])
    output:
        protected("results/called/{sample}_raw_snps_indels.vcf")
    params:
        tdir = config['TEMPDIR'],
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
        # partition = config['PARTITION']['CPU']
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bams}"
    shell:
        'gatk HaplotypeCaller '
        '--java-options "-Xmx{resources.mem_gb}G" '
        '-I {input.bams} '
        '-R {input.refgenome} '
        # '-D {input.dbsnp} '
        '-O {output} '
        '--tmp-dir {params.tdir} {params.padding} {params.intervals} '
        '&> {log}'


rule gatk_BaseRecalibrator:
    input:
        bams = "results/mapped/{sample}_sorted_mkdups.bam",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
    output:
        report("results/mapped/{sample}_recalibration_report.grp", caption = "report/recalibration.rst", category = "Base recalibration")
    params:
        tdir = config['TEMPDIR'],
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
        # partition = config['PARTITION']['CPU']
    message:
        "Generating a recalibration table for {input.bams}"
    shell:
        'gatk BaseRecalibrator '
        '--java-options "-Xmx{resources.mem_gb}G" '
        '-I {input.bams} '
        '-R {input.refgenome} '
        '-O {output} '
        '--tmp-dir {params.tdir} {params.padding} {params.intervals} {params.recalibration_resources} '
        '&> {log}'


rule gatk_ApplyBQSR:
    input:
        bam = "results/mapped/{sample}_sorted_mkdups.bam",
        recal = "results/mapped/{sample}_recalibration_report.grp",
        refgenome = expand("{refgenome}", refgenome = config['REFGENOME'])
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
        # partition = config['PARTITION']['CPU']
    message:
        "Applying base quality score recalibration and producing a recalibrated BAM file for {input.bam}"
    shell:
        'gatk ApplyBQSR '
        '--java-options "-Xmx{resources.mem_gb}G" '
        '-I {input.bam} '
        '-bqsr {input.recal} '
        '-R {input.refgenome} '
        '-O {output} {params.padding} {params.intervals}'


