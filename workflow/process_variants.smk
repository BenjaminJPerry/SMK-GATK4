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
        "results/02_snvs/merged.filteredsnvs.MQM60.freebayes.vcf.gz",
        "results/02_snvs/merged.filteredsnvs.MQM60.bcftools.vcf.gz",


rule freebayes_vcf: #ADDED TEMPORARILY TODO REMOVE AGAIN
    priority: 1000
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        vcf = temp("results/02_snvs/{samples}.rawsnvs.freebayes.vcf"),
    log:
        "logs/freebayes_vcf.{samples}.log"
    benchmark:
        "benchmarks/freebayes_vcf.{samples}.tsv"
    threads: 2
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "freebayes "
        "--standard-filters "
        "--pooled-continuous "
        "--trim-complex-tail "
        "-F 0.01 "
        "-f {input.referenceGenome} {input.bam} > {output.vcf}"


rule bgzip_freebayes_vcf:
    priority:100
    input:
        vcf = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf",
    output:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz",
    benchmark:
        "benchmarks/bgzip_freebayes_vcf.{samples}.tsv"
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
        
        bgzip -c -l 8 --threads {threads} {input.vcf} > {output.vcfgz}

        """


rule index_freebayes_vcf:
    priority:100
    input:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz",
    output:
        csi = temp("results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz.csi"),
    benchmark:
        "benchmarks/index_freebayes_vcf.{samples}.tsv"
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
        
        bcftools index --threads {threads} {input.vcfgz} -o {output.csi}

        """


rule index_bcftools_vcf:
    priority:100
    input:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz",
    output:
        csi = temp("results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz.csi"),
    benchmark:
        "benchmarks/index_bcftools_vcf.{samples}.tsv"
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
        
        bcftools index --threads {threads} {input.vcfgz} -o {output.csi}

        """


rule merge_freebayes_vcf: #TODO
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz"
    benchmark:
        "benchmarks/merge_freebayes_vcf.tsv"
    threads: 16
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools merge --threads {threads} {input.vcfgz} -Oz8 -o {output.merged}

        """


rule merge_bcftools_vcf: #TODO
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz"
    benchmark:
        "benchmarks/merge_bcftools_vcf.tsv"
    threads: 16
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools merge --threads {threads} {input.vcfgz} -Oz8 -o {output.merged}

        """


rule index_freebayes_merged:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz"
    output:
        index = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz.tbi"
    benchmark:
        "benchmarks/filter_freebayes_vcf.tsv"
    threads: 16
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "tabix -f results/02_snvs/merged.rawsnvs.freebayes.vcf.gz"


rule index_bcftools_merged:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz"
    output:
        index = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz.tbi"
    benchmark:
        "benchmarks/index_bcftools_merged.tsv"
    threads: 16
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "tabix -f {input.merged} "


rule filter_freebayes_vcf:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz",
        index = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz.tbi"
    output:
        filtered = "results/02_snvs/merged.filteredsnvs.MQM60.freebayes.vcf.gz"
    benchmark:
        "benchmarks/filter_freebayes_vcf.tsv"
    threads: 16
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "vcffilter -f 'MQM = 60' {input.merged} > {output.filtered} "


rule filter_bcftools_vcf:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz",
        index = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz.tbi"
    output:
        filtered = "results/02_snvs/merged.filteredsnvs.MQM60.bcftools.vcf.gz"
    benchmark:
        "benchmarks/filter_bcftools_vcf.tsv"
    threads: 16
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "vcffilter -f 'MQM = 60' {input.merged} > {output.filtered} "

