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
CHROM = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28', 'chr29') # , 'chrY', 'chrX'


# wildcard_constraints:
#    chromosome = '\w+',
#    samples = '\w+'


rule all:
    input:
        "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz",
        "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz.csi"



rule gatk_HaplotypeCaller_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        vcf_chrom = temp("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz.csi")
    params:
        chromosome = '{chromosome}'
    log:
        "logs/gatk_HaplotypeCaller_vcf.{samples}.{chromosome}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_vcf.{samples}.{chromosome}.tsv"
    conda:
        "bcftools"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 12),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 240),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.4.0.0-gimkl-2022a; '
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
        '&> {log}.attempt.{resources.attempt} && '
        'rm results/02_snvs/{wildcards.samples}.rawsnvs.{wildcards.chromosome}.haplotypeCaller.vcf.gz.tbi; '
        'bcftools index --threads {threads} {output.vcf_chrom} -o {output.csi}'


rule concatenate_replicons_vcf: #TODO
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz", chromosome = CHROM, allow_missing=True),
        csi = expand("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz", chromosome = CHROM, allow_missing=True),
    output:
        merged = temp("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi"),
    benchmark:
        "benchmarks/{samples}_concatenate_replicons_vcf.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools concat --threads {threads} {input.vcfgz} | bcftools sort -O z8 -o {output.merged} - ;

        bcftools index --threads {threads} {output.merged} -o {output.csi}

        """


rule merge_animals_vcf:
    priority:100
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz", samples = SAMPLES),
        csi = expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.rawsnvs.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/merged.rawsnvs.haplotypeCaller.vcf.gz.csi",
    benchmark:
        "benchmarks/merge_animals_vcf.tsv"
    threads: 16
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools merge --threads {threads} {input.vcfgz} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 


        """


rule view_haplotype_chrom:
    priority:100
    input:
        merged_vcf = "results/02_snvs/merged.rawsnvs.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/merged.rawsnvs.haplotypeCaller.vcf.gz.csi",
    output:
        filtered_vcf = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz",
        filtered_vcf_csi = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz.csi"
    params:
        chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29" #TODO move to config; Also, removed ChrM ,chrY, chrX
    benchmark:
        "benchmarks/view_haplotypeCaller_chrom.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.merged_vcf} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes};

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi}; 

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """

