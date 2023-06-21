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


SAMPLES, = glob_wildcards("results/01_mapping/{samples}.sorted.mkdups.merged.bam")

AHJNKHDSX5,= glob_wildcards("fastq/AHJNKHDSX5/{samples}_R1.fastq.gz")
BHJNVTDSX5, = glob_wildcards("fastq/BHJNVTDSX5/{samples}_R1.fastq.gz")
MGI, = glob_wildcards("fastq/MGI/{samples}_R1.fastq.gz")

rule all:
    input:
        "results/00_stats/seqkit.fastq.stats.txt",
        expand("results/00_stats/fastqc/{samples}_R2_fastqc.zip", samples = AHJNKHDSX5),
        expand("results/00_stats/fastqc/{samples}_R2_fastqc.zip", samples = BHJNVTDSX5),
        expand("results/00_stats/fastqc/{samples}_R2_fastqc.zip", samples = MGI),
        expand("results/00_stats/fastqc/{samples}_R1_fastqc.zip", samples = AHJNKHDSX5),
        expand("results/00_stats/fastqc/{samples}_R1_fastqc.zip", samples = BHJNVTDSX5),
        expand("results/00_stats/fastqc/{samples}_R1_fastqc.zip", samples = MGI),
        expand("results/00_stats/{samples}.sorted.mkdups.merged.bam.samtools_stats.txt", samples = SAMPLES),
        expand("results/00_stats/{samples}.mosdepth.summary.txt", samples = SAMPLES),


rule fastqc_AHJNKHDSX5:
    priority: 100
    input:
        fastq1 = "fastq/AHJNKHDSX5/{samples}_R1.fastq.gz",
        fastq2 = "fastq/AHJNKHDSX5/{samples}_R2.fastq.gz",
    output:
        zip1 = "results/00_stats/fastqc/{samples}_R1_fastqc.zip",
        zip2 = "results/00_stats/fastqc/{samples}_R2_fastqc.zip"
    conda:
        "fastqc"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,    
    shell:
        "fastqc "
        "-o results/00_stats/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.fastq1} {input.fastq2}"


rule fastqc_BHJNVTDSX5:
    priority: 100
    input:
        fastq1 = "fastq/BHJNVTDSX5/{samples}_R1.fastq.gz",
        fastq2 = "fastq/BHJNVTDSX5/{samples}_R2.fastq.gz",
    output:
        zip1 = "results/00_stats/fastqc/{samples}_R1_fastqc.zip",
        zip2 = "results/00_stats/fastqc/{samples}_R2_fastqc.zip"
    conda:
        "fastqc"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,    
    shell:
        "fastqc "
        "-o results/00_stats/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.fastq1} {input.fastq2}"


rule fastqc_MGI:
    priority: 100
    input:
        fastq1 = "fastq/MGI/{samples}_R1.fastq.gz",
        fastq2 = "fastq/MGI/{samples}_R2.fastq.gz",
    output:
        zip1 = "results/00_stats/fastqc/{samples}_R1_fastqc.zip",
        zip2 = "results/00_stats/fastqc/{samples}_R2_fastqc.zip"
    conda:
        "fastqc"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,    
    shell:
        "fastqc "
        "-o results/00_stats/fastqc/ "
        "-q "
        "-t {threads} "
        "{input.fastq1} {input.fastq2}"



rule seqkit_stats:
    priority: 1000
    input:
        fastqA1 = expand("fastq/AHJNKHDSX5/{samples}_R1.fastq.gz", samples = AHJNKHDSX5),
        fastqA2 = expand("fastq/AHJNKHDSX5/{samples}_R2.fastq.gz", samples = AHJNKHDSX5),
        fastqB1 = expand("fastq/BHJNVTDSX5/{samples}_R1.fastq.gz", samples = BHJNVTDSX5),
        fastqB2 = expand("fastq/BHJNVTDSX5/{samples}_R2.fastq.gz", samples = BHJNVTDSX5),
        fastqM1 = expand("fastq/MGI/{samples}_R1.fastq.gz", samples = MGI),
        fastqM2 = expand("fastq/MGI/{samples}_R2.fastq.gz", samples = MGI),
    output:
        "results/00_stats/seqkit.fastq.stats.txt"
    benchmark:
        "benchmarks/seqkit_stats.txt"
    conda:
        "seqkit"
    threads: 64
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 360 + ((attempt - 1) * 600),
        partition="large,milan",
    shell:
        "seqkit stats -j {threads} -a {input.fastqA1} {input.fastqA2} {input.fastqB1} {input.fastqB2} {input.fastqM1} {input.fastqM2} > {output} "



rule samtools_stats_merged:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa"
    output:
        stats = "results/00_stats/{samples}.sorted.mkdups.merged.bam.samtools_stats.txt"
    log:
        "logs/samtools_stats_merged.{samples}.log"
    benchmark:
        "benchmarks/samtools_stats_merged.{samples}.tsv"
    conda:
        "bwa"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "samtools stats --threads {threads} -r {input.referenceGenome} {input.bam} > {output.stats} "


rule mosdepth_stats_merged:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
    output:
        stats = "results/00_stats/{samples}.mosdepth.summary.txt",
    log:
        "logs/mosdepth_stats_merged.{samples}.log"
    benchmark:
        "benchmarks/mosdepth_stats_merged.{samples}.tsv"
    conda:
        "mosdepth"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "mosdepth "
        "--fast-mode " # dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases).
        "--no-per-base " # dont output per-base depth.
        "--threads {threads} "
        "results/00_stats/{wildcards.samples} " # output prefix
        "{input.bam} "


# rule bcftools_stats_freebayes_raw: #TODO
#     priority: 100
#     input:
#         bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
#         referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
#     output:
#         vcf = temp("results/02_snvs/{samples}.rawsnvs.freebayes.vcf"),
#     log:
#         "logs/freebayes_vcf.{samples}.log"
#     benchmark:
#         "benchmarks/freebayes_vcf.{samples}.tsv"
#     threads: 2
#     conda:
#         "freebayes"
#     resources:
#         mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
#         time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
#         partition = "large,milan",
#         DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
#         attempt = lambda wildcards, attempt: attempt,
#     shell:
#         "freebayes "
#         "--standard-filters "
#         "--pooled-continuous "
#         "--trim-complex-tail "
#         "-F 0.01 "
#         "-f {input.referenceGenome} {input.bam} > {output.vcf}"

# rule bcftools_stats_bcftools_raw:
# rule bcftools_stats_freebayes_filtered:
# rule bcftools_stats_bcftools_filtered:

# rule rtg_vcfstats_freebayes_raw:
# rule rtg_vcfstats_freebayes_filtered:
# rule rtg_vcfstats_bcftools_raw:
# rule rtg_vcfstats_bcftools_filtered:


