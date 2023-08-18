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
        "results/02_snvs/merged.chrom.DPFilt.bcftools.vcf.gz",
        "results/02_snvs/merged.chrom.DPFilt.freebayes.vcf.gz",

        # "results/02_snvs/merged.chrom.bcftools.vcf.gz.csi",
        # "results/02_snvs/merged.chrom.freebayes.vcf.gz.csi"
        
        #"results/02_snvs/merged.rawsnvs.bcftools.vcf.gz",
        #"results/02_snvs/merged.rawsnvs.freebayes.vcf.gz",

        #expand("results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz", samples = SAMPLES),
        #expand("results/02_snvs/{samples}.rawsnvs.varscan2.vcf.gz", samples = SAMPLES),
        #expand("results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz", samples = SAMPLES),
        #expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz", samples = SAMPLES),


### bcftools 


rule bcftools_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        vcf = temp("results/02_snvs/{samples}.rawsnvs.DPFilt.bcftools.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.DPFilt.bcftools.vcf.gz.csi"),
    log:
        "logs/bcftools_vcf.{samples}.log"
    benchmark:
        "benchmarks/bcftools_vcf.{samples}.tsv"
    threads: 24
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools mpileup --seed 1953 --threads {threads} "
        "--max-depth 250 " # Max raw per-file depth; avoids excessive memory usage [250]
        "-q 30 " # skip alignment with mapQ less than
        "-Q 20 " # Skip bases with baseQ/BAQ less than
        "-m 10 " # Minimum number gapped reads for indel candidates
        "-f {input.referenceGenome} "
        "{input.bam} "
        "| bcftools call "
        "--threads {threads} "
        "-v " # Output variant sites only
        "-m " # Alternative model for multiallelic and rare-variant calling
        "| bcftools view --threads {threads} -O z8 -e 'INFO/DP<10 || INFO/DP>500' -o {output.vcf}; "
        "bcftools index --threads {threads} {output.vcf} -o  {output.csi} "


rule merge_bcftools_vcf:
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.DPFilt.bcftools.vcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.DPFilt.bcftools.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.rawsnvs.DPFilt.bcftools.vcf.gz",
        csi = "results/02_snvs/merged.rawsnvs.DPFilt.bcftools.vcf.gz.csi"
    benchmark:
        "benchmarks/merge_bcftools_vcf.tsv"
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
        
        bcftools merge --threads {threads} {input.vcfgz} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 


        """


### freebayes


rule freebayes_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        vcfgz = temp("results/02_snvs/{samples}.rawsnvs.DPFilt.freebayes.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.DPFilt.freebayes.vcf.gz.csi"),
    log:
        "logs/freebayes_vcf.{samples}.log"
    benchmark:
        "benchmarks/freebayes_vcf.{samples}.tsv"
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
        "--limit-coverage 250 " # Match the other variant callers
        #"--pooled-continuous " # Output all alleles which pass input filters, regardles of genotyping outcome or model.
        #"--trim-complex-tail " # Trim complex tails.
        #"-F 0.01 " # minimum fraction of observations supporting alternate allele within one individual [0.05]
        "-f {input.referenceGenome} {input.bam} "
        "| bcftools view --threads {threads} -O z8 -e 'INFO/DP<10 || INFO/DP>500' -o {output.vcfgz}; "
        "bcftools index --threads {threads} {output.vcfgz} -o  {output.csi} "


rule merge_freebayes_vcf:
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.DPFilt.freebayes.vcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.DPFilt.freebayes.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.rawsnvs.DPFilt.freebayes.vcf.gz",
        csi = "results/02_snvs/merged.rawsnvs.DPFilt.freebayes.vcf.gz.csi",
    benchmark:
        "benchmarks/merge_freebayes_vcf.tsv"
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
        
        bcftools merge --threads {threads} {input.vcfgz} -O z8 -o {output.merged} &&

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 

        """


rule view_bcftools_chrom:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.DPFilt.bcftools.vcf.gz",
        csi = "results/02_snvs/merged.rawsnvs.DPFilt.bcftools.vcf.gz.csi"
    output:
        filtered_vcf = "results/02_snvs/merged.chrom.DPFilt.bcftools.vcf.gz",
        filtered_vcf_csi = "results/02_snvs/merged.chrom.DPFilt.bcftools.vcf.gz.csi"
    params:
        chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chrX,chrY" #TODO move to config; Also, removed ChrM
    benchmark:
        "benchmarks/view_bcftools_chrom.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 720 + ((attempt - 1) * 720),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.merged_vcf} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes} &&

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi} &&

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """


rule view_freebayes_chrom:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.DPFilt.freebayes.vcf.gz",
        csi = "results/02_snvs/merged.rawsnvs.DPFilt.freebayes.vcf.gz.csi",
    output:
        filtered_vcf = "results/02_snvs/merged.chrom.DPFilt.freebayes.vcf.gz",
        filtered_vcf_csi = "results/02_snvs/merged.chrom.DPFilt.freebayes.vcf.gz.csi"
    params:
        chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chrX,chrY" #TODO move to config; Also, removed ChrM
    benchmark:
        "benchmarks/view_freebayes_chrom.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 720 + ((attempt - 1) * 720),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.merged_vcf} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes} &&

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi} && 

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """

