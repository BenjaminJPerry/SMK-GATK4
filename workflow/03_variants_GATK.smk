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
CHROM = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28', 'chr29', 'chrX', 'chrY', 'chrM')


# wildcard_constraints:
#    chromosome = '\w+',
#    samples = '\w+'


rule all:
    input:
        "results/02_snvs/merged.rawsnvs.haplotypeCaller.vcf.gz",
        "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz.csi"



rule gatk_HaplotypeCaller_vcf:
    priority: 1
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
        #chromosome = '{chromosome}'
    output:
        vcf_chrom = temp("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz"),
        vcf_chrom_index = temp("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz.tbi")
    params:
        chromosome = '{chromosome}'
    log:
        "logs/gatk_HaplotypeCaller_vcf.{samples}.{chromosome}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_vcf.{samples}.{chromosome}.tsv"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 720 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.4.0.0-gimkl-2022a ; '
        'gatk --java-options "-Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}"  '
        'HaplotypeCaller '
        '--base-quality-score-threshold 20 ' 
        '--min-base-quality-score 20 '
        '--create-output-variant-index '
        '-R {input.referenceGenome} '
        '-I {input.bam} '
        '-L {params.chromosome} '
        '-O {output.vcf_chrom} '
        '--tmp-dir {resources.DTMP} '
        '&> {log}.attempt.{resources.attempt} '


rule index_replicons_vcf:
    priority:100
    input:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz",
        vcf_chrom_index = "results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz.tbi"
    output:
        csi = temp("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz.csi"),
    benchmark:
        "benchmarks/index_bcftools_vcf.{samples}.{chromosome}.tsv"
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


rule concatenate_replicons_vcf: #TODO
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz", chromosome = CHROM, allow_missing=True),
        csi = expand("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz.csi", chromosome = CHROM, allow_missing=True),
    output:
        merged = temp("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz")
    benchmark:
        "benchmarks/{samples}_merge_bcftools_vcf.tsv"
    threads: 8
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
        
        bcftools concat --threads {threads} {input.vcfgz} | bcftools sort -Oz8 -o {output.merged} - 

        """



rule index_animals_vcf:
    priority:100
    input:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz",
    output:
        csi = temp("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi"),
    benchmark:
        "benchmarks/index_varscan2_vcf.{samples}.tsv"
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


rule merge_animals_vcf: #TODO
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz.csi",
    benchmark:
        "benchmarks/merge_varscan2_vcf.tsv"
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
        
        bcftools merge --threads {threads} {input.vcfgz} -Oz8 -o {output.merged} &&

        bcftools index --threads {threads} {output.merged} -o {output.csi} &&

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 


        """

