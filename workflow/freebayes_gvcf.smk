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
        "results/02_snvs/merged.QUAL20.freebayes.gvcf.gz.pigmentSNPs.vcf",


rule freebayes_gvcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        gvcf = "results/02_snvs/{samples}.QUAL20.freebayes.gvcf",
    log:
        "logs/freebayes_gvcf.{samples}.log"
    benchmark:
        "benchmarks/freebayes_gvcf.{samples}.tsv"
    threads: 2
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "freebayes "
        "--standard-filters " # Use stringent input base and mapping quality filters. Equivalent to -m 30 -q 20 -R 0 -S 0
        "--pooled-continuous " # Output all alleles which pass input filters, regardles of genotyping outcome or model.
        "--trim-complex-tail " # Trim complex tails.
        "--gvcf "
        "-F 0.01 " # minimum fraction of observations supporting alternate allele within one individual [0.05]
        "-f {input.referenceGenome} {input.bam} | "
        "vcffilter -f 'QUAL > 20' "
        "> {output.gvcf}  "


rule bgzip_freebayes_gvcf:
    priority:100
    input:
        gvcf = "results/02_snvs/{samples}.QUAL20.freebayes.gvcf",
    output:
        gvcfgz = "results/02_snvs/{samples}.QUAL20.freebayes.gvcf.gz",
    benchmark:
        "benchmarks/bgzip_freebayes_gvcf.{samples}.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 16 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 240),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bgzip -c -l 8 --threads {threads} {input.gvcf} > {output.gvcfgz}

        """


rule index_freebayes_gvcf:
    priority:100
    input:
        gvcfgz = "results/02_snvs/{samples}.QUAL20.freebayes.gvcf.gz",
    output:
        csi = temp("results/02_snvs/{samples}.QUAL20.freebayes.gvcf.gz.csi"),
    benchmark:
        "benchmarks/index_freebayes_gvcf.{samples}.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 16 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 240),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools index --threads {threads} {input.vcfgz} -o {output.csi}

        """


rule merge_freebayes_gvcf: #TODO
    priority:100
    input:
        gvcfgz = expand("results/02_snvs/{samples}.QUAL20.freebayes.gvcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.freebayes.gvcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.rawsnvs.freebayes.gvcf.gz"
    benchmark:
        "benchmarks/merge_freebayes_gvcf.tsv"
    threads: 16
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools merge --threads {threads} {input.gvcfgz} -Oz8 -o {output.merged}

        """


rule bcftools_view_freebayes_regions:
    priority:100
    input:
        gvcf = "results/02_snvs/merged.QUAL20.freebayes.gvcf.gz",
        regions = "resources/snp_targets.txt"
    output:
        filtered_snps = "results/02_snvs/merged.QUAL20.freebayes.gvcf.gz.pigmentSNPs.vcf"
    threads: 2
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 6 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        " bcftools view "
        "-R {input.regions} "
        "--threads {threads} "
        "{input.filtered} "
        "> {output.filtered_snps} "
