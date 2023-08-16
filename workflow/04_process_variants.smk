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
        "results/03_filtered/merged.chrom.freebayes.QUAL60.vcf.gz",
        "results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz",
        "results/03_filtered/merged.chrom.haplotypeCaller.QUAL60.vcf.gz",

        #"results/03_filtered/merged.chrom.freebayes.QUAL60.vcf.gz.pigmentSNPs.vcf",
        #"results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz.pigmentSNPs.vcf",
        #"results/03_filtered/merged.chrom.haplotypeCaller.QUAL60.vcf.gz.pigmentSNPs.vcf",



rule filter_freebayes_vcf_QUAL60:
    priority:100
    input:
        merged = "results/02_snvs/merged.chrom.freebayes.vcf.gz",
        index = "results/02_snvs/merged.chrom.freebayes.vcf.gz.csi"
    output:
        filtered = "results/03_filtered/merged.chrom.freebayes.QUAL60.vcf.gz",
        csi = "results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz.csi"
    benchmark:
        "benchmarks/filter_freebayes_vcf_QUAL60.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''
        sleep 5
        bcftools view -e 'QUAL<60' {input.merged} -O z8 -o {output.filtered} &&
        bcftools index --threads {threads} {output.filtered} -o {output.csi} 

        echo "Total snps in {input.merged} at QUAL>=60: $(bcftools view --threads {threads} -i 'QUAL>=60' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=50: $(bcftools view --threads {threads} -i 'QUAL>=50' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=40: $(bcftools view --threads {threads} -i 'QUAL>=40' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=30: $(bcftools view --threads {threads} -i 'QUAL>=30' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=20: $(bcftools view --threads {threads} -i 'QUAL>=20' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt &&

        exit 0;
        '''

rule filter_bcftools_vcf_QUAL60:
    priority:100
    input:
        merged = "results/02_snvs/merged.chrom.bcftools.vcf.gz",
        index = "results/02_snvs/merged.chrom.bcftools.vcf.gz.csi"
    output:
        filtered = "results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz",
        csi = "results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz.csi"
    benchmark:
        "benchmarks/filter_bcftools_vcf_QUAL60.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''
        sleep 5
        bcftools view -e 'QUAL<60' {input.merged} -O z8 -o {output.filtered} &&
        bcftools index --threads {threads} {output.filtered} -o {output.csi} 

        echo "Total snps in {input.merged} at QUAL>=60: $(bcftools view --threads {threads} -i 'QUAL>=60' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=50: $(bcftools view --threads {threads} -i 'QUAL>=50' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=40: $(bcftools view --threads {threads} -i 'QUAL>=40' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=30: $(bcftools view --threads {threads} -i 'QUAL>=30' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {input.merged} at QUAL>=20: $(bcftools view --threads {threads} -i 'QUAL>=20' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt &&

        exit 0;

        '''


rule filter_haplotypeCaller_vcf_QUAL60:
    priority:100
    input:
        merged = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz",
        index = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz.csi"
    output:
        filtered = "results/03_filtered/merged.chrom.haplotypeCaller.QUAL60.vcf.gz",
        csi = "results/03_filtered/merged.chrom.haplotypeCaller.QUAL60.vcf.gz.csi"
    benchmark:
        "benchmarks/filter_haplotypeCaller_vcf_QUAL60.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''
        sleep 5
        bcftools view -e 'QUAL<60' {input.merged} -O z8 -o {output.filtered} &&
        bcftools index --threads {threads} {output.filtered} -o {output.csi} 

        echo "Total snps in {input.merged} at QUAL>=60: $(bcftools view --threads {threads} -i 'QUAL>=60' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 
        echo "Total snps in {input.merged} at QUAL>=50: $(bcftools view --threads {threads} -i 'QUAL>=50' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;
        echo "Total snps in {input.merged} at QUAL>=40: $(bcftools view --threads {threads} -i 'QUAL>=40' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 
        echo "Total snps in {input.merged} at QUAL>=30: $(bcftools view --threads {threads} -i 'QUAL>=30' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 
        echo "Total snps in {input.merged} at QUAL>=20: $(bcftools view --threads {threads} -i 'QUAL>=20' {input.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt &&

        exit 0;

        '''



rule view_bcftools_regions: #TODO Updated files
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz",
        regions = "resources/snp_targets_100bp.txt"
    output:
        filtered_snps = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz.pigmentSNPs.vcf"
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


rule view_freebayes_regions: #TODO Updated files
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf.gz",
        regions = "resources/snp_targets_100bp.txt"
    output:
        filtered_snps = "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf.gz.pigmentSNPs.vcf"
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


rule view_haplotypeCaller_regions: #TODO Updated files
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.haplotypeCaller.vcf.gz",
        regions = "resources/snp_targets_100bp.txt"
    output:
        filtered_snps = "results/03_filtered/merged.filteredsnvs.QUAL60.haplotypeCaller.vcf.gz.pigmentSNPs.vcf"
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



rule isec_bcftools_LIC565:
    priority:100
    input:
        filtered = "results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz",
        csi = "results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz.csi",
        snps_LIC = "resources/LIC_565.ch.frmt.sorted.vcf.gz",
        snps_LIC_csi = "resources/LIC_565.ch.frmt.sorted.vcf.gz.csi",
    output:
        LICFiltered = "results/03_filtered/merged.chrom.bcftools.QUAL60.LIC565.vcf.gz",
        csi = "results/03_filtered/merged.chrom.bcftools.QUAL60.LIC565.vcf.gz.csi"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/03_filtered/isec_bcftools_LIC565 {input.snps_LIC} {input.filtered} &&

        #TODO view to vcf.gz

        rm -r results/03_filtered/isec_bcftools_LIC565 &&

        bcftools index --threads {threads} {output.LICFiltered} -o {output.csi} && 

        echo "Total snps in {output.LICFiltered}: $(bcftools view --threads 6 {output.LICFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """

rule isec_freebayes_LIC565:


rule isec_haplotypeCaller_LIC565:



rule isec_bcftools_TBulls:


rule isec_freebayes_TBulls:


rule isec_haplotypeCaller_TBulls:



rule filter_DP_bcftools:


rule filter_DP_freebayes:


rule filter_DP_haplotypeCaller:



