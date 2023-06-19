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


rule all:
    input:
        expand("results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz", samples = SAMPLES),
        expand("results/02_snvs/{samples}.rawsnvs.freebayes.vcf", samples = SAMPLES),
        expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz", samples = SAMPLES),
        expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.gvcf.gz", samples = SAMPLES),


rule bgzip_freebayes_vcf: #TODO
    priority:100
    input:
        vcf = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf",
    output:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz",
    benchmark:
        "benchmarks/bgzip_freebayes_vcf"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 16 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 240),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bgzip -c -i -l 8 --threads {threads} {input.vcf} > {output.vcfgz}

        """


rule index_bcftools_vcf: #TODO
    priority:100
    input:
        vcf = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz",
    output:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz.csi",
    benchmark:
        "benchmarks/bgzip_freebayes_vcf"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 16 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 240),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bgzip -c -i -l 8 --threads {threads} {input.vcf} > {output.vcfgz}

        """



rule merge_bcftools_vcf: #TODO
    priority:100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        gvcf = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.gvcf.gz",
    log:
        "logs/gatk_HaplotypeCaller.gvcf.{samples}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller.gvcf.{samples}.tsv"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 7200 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.3.0.0-gimkl-2022a ; '
        'gatk --java-options "-Xms2G -Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'HaplotypeCaller '
        '--min-base-quality-score 30 '
        '-I {input.bam} '
        '-R {input.referenceGenome} '
        '-O {output.gvcf} '
        '-ERC GVCF '
        '--tmp-dir {resources.DTMP} '
        '&> {log}.attempt.{resources.attempt} '


rule merge_freebayes_vcf: #TODO
    priority:100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        gvcf = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.gvcf.gz",
    log:
        "logs/gatk_HaplotypeCaller.gvcf.{samples}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller.gvcf.{samples}.tsv"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 7200 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.3.0.0-gimkl-2022a ; '
        'gatk --java-options "-Xms2G -Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'HaplotypeCaller '
        '--min-base-quality-score 30 '
        '-I {input.bam} '
        '-R {input.referenceGenome} '
        '-O {output.gvcf} '
        '-ERC GVCF '
        '--tmp-dir {resources.DTMP} '
        '&> {log}.attempt.{resources.attempt} '