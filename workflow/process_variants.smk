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
        "results/03_filtered/merged.filteredsnvs.QUAL20.freebayes.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL20.bcftools.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL30.freebayes.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL30.bcftools.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL40.freebayes.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL40.bcftools.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL50.freebayes.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL50.bcftools.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf.gz",
        "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz",
        "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf.gz.pigmentSNPs.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz.pigmentSNPs.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL20.freebayes.vcf.gz.pigmentSNPs.vcf",
        "results/03_filtered/merged.filteredsnvs.QUAL20.bcftools.vcf.gz.pigmentSNPs.vcf"


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
        partition = "large,milan",
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
        partition = "large,milan",
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
        partition = "large,milan",
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
        partition = "large,milan",
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
        partition = "large,milan",
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
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "tabix -f {input.merged} "


rule filter_freebayes_vcf_QUAL60:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz",
        index = "results/02_snvs/merged.rawsnvs.freebayes.vcf.gz.tbi"
    output:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf"
    benchmark:
        "benchmarks/filter_freebayes_vcf.tsv"
    threads: 16
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "vcffilter -f 'QUAL > 60' {input.merged} > {output.filtered} " #TODO COMPRESS OUTPUT


rule filter_bcftools_vcf_QUAL60:
    priority:100
    input:
        merged = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz",
        index = "results/02_snvs/merged.rawsnvs.bcftools.vcf.gz.tbi"
    output:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf"
    benchmark:
        "benchmarks/filter_bcftools_vcf.tsv"
    threads: 16
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "vcffilter -f 'QUAL > 60' {input.merged} > {output.filtered} " #TODO COMPRESS OUTPUT


rule bcftools_view_freebayes_fvcf: #TODO
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf",
        regions = "resources/regions.txt"
    output:
        regionSNVs = "results/02_snvs/merged.rawsnvs.bcftools.vcf"
    benchmark:
        "benchmarks/bcftools_view_freebayes_fvcf.tsv"
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
        " bcftools view "
        "--threads {threads} "
        "-R <regions.txt> "
        "--with-header "
        "{input.vcfgz} "
        "-Oz8 -o {output.merged} "


rule bcftools_view_bcf_freebayes:
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf",
    output:
        bcf = "results/03_filtered/merged.filteredsnvs.QUAL60.freebayes.vcf.gz"
    threads: 32
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        " bcftools view "
        "--threads {threads} "
        "{input.filtered} "
        "-Oz8 -o {output.bcf} && "
        "bcftools index --threads {threads} {output.bcf} "


rule bcftools_view_bcf_bcftools:
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf",
    output:
        bcf = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz"
    threads: 32
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        " bcftools view "
        "--threads {threads} "
        "{input.filtered} "
        "-Oz8 -o {output.bcf} && "
        "bcftools index --threads {threads} {output.bcf} "


rule bcftools_view_bcftools_regions:
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


rule bcftools_view_freebayes_regions:
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


rule bcftools_view_bcf_freebayes_QUAL20:
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL20.freebayes.vcf",
    output:
        bcf = "results/03_filtered/merged.filteredsnvs.QUAL20.freebayes.vcf.gz"
    threads: 32
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        " bcftools view "
        "--threads {threads} "
        "{input.filtered} "
        "-Oz8 -o {output.bcf} && "
        "bcftools index --threads {threads} {output.bcf} "


rule bcftools_view_bcf_bcftools_QUAL20:
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL20.bcftools.vcf",
    output:
        bcf = "results/03_filtered/merged.filteredsnvs.QUAL20.bcftools.vcf.gz"
    threads: 32
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 1440),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        " bcftools view "
        "--threads {threads} "
        "{input.filtered} "
        "-Oz8 -o {output.bcf} && "
        "bcftools index --threads {threads} {output.bcf} "



rule bcftools_view_bcftools_regions_QUAL20:
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL20.bcftools.vcf.gz",
        regions = "resources/snp_targets_100bp.txt"
    output:
        filtered_snps = "results/03_filtered/merged.filteredsnvs.QUAL20.bcftools.vcf.gz.pigmentSNPs.vcf"
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


rule bcftools_view_freebayes_regions_QUAL20:
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL20.freebayes.vcf.gz",
        regions = "resources/snp_targets_100bp.txt"
    output:
        filtered_snps = "results/03_filtered/merged.filteredsnvs.QUAL20.freebayes.vcf.gz.pigmentSNPs.vcf"
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


rule vcftools_vcfstats_bcftools_fvcf: #TODO


rule rtg_vcfstats__bcftools_fvcf: #TODO


rule bcftools_stats_bcftools_fvcf: #TODO