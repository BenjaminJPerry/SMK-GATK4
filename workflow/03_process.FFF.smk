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


# ENTRY POINTS
        # "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.vcf.gz",
        # "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.vcf.gz",
        # "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz",


rule all:
    input:
        "results/05_ensemble/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.intersect.vcf.gz",
        "results/05_ensemble/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.eva.intersect.vcf.gz",
        "results/05_ensemble/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.eva.intersect.vcf.gz",


rule get_eva_snvs:
    output:
        vcf = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz",
        csi = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz.csi",
    params:
        vcf = "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/ovis_aries/ARSUI_Ramb_v2.0/9940_GCA_016772045.1_current_ids.vcf.gz",
        index = "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/ovis_aries/ARSUI_Ramb_v2.0/9940_GCA_016772045.1_current_ids.vcf.gz.csi",
    threads: 2
    shell:
        """
        
        wget -t 5 -O {output.vcf} {params.vcf} ;
        wget -t 5 -O {output.csi} {params.index} ; 
        
        """


rule rename_eva_snvs:
    input:
        vcf = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz",
        csi = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz.csi",
        sed_file = "resources/eva_rename.sed",
    output:
        vcf = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi",
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        zcat {input.vcf} | sed -f {input.sed_file} | bcftools view --threads {threads} -O z8 -o {output.vcf} -

        bcftools index --threads {threads} {output.vcf} -o {output.csi} ; 
        
        """


rule isec_bcftools_eva:
    priority:100
    input:
        filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi",
        eva = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi"
    output:
        eva_filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz.csi",
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/04_merged/isec_bcftools_eva {input.eva} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/04_merged/isec_bcftools_eva/0001.vcf -o {output.eva_filtered} ;

        bcftools index --threads {threads} {output.eva_filtered} -o {output.csi} ; 

        rm -r results/04_merged/isec_bcftools_eva ;

        echo "Total snps in {output.eva_filtered}: $(bcftools view --threads {threads} {output.eva_filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule isec_freebayes_eva:
    priority:100
    input:
        filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.vcf.gz.csi",
        eva = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi"
    output:
        eva_filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.eva.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.eva.vcf.gz.csi",
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/04_merged/isec_freebayes_eva {input.eva} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/04_merged/isec_freebayes_eva/0001.vcf -o {output.eva_filtered} ;

        bcftools index --threads {threads} {output.eva_filtered} -o {output.csi} ; 

        rm -r results/04_merged/isec_freebayes_eva ;

        echo "Total snps in {output.eva_filtered}: $(bcftools view --threads {threads} {output.eva_filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule isec_haplotypeCaller_eva:
    priority:100
    input:
        filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi",
        eva = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi"
    output:
        eva_filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.eva.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.eva.vcf.gz.csi",
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/04_merged/isec_haplotypeCaller_eva {input.eva} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/04_merged/isec_haplotypeCaller_eva/0001.vcf -o {output.eva_filtered} ;

        bcftools index --threads {threads} {output.eva_filtered} -o {output.csi} ; 

        rm -r results/04_merged/isec_haplotypeCaller_eva ;

        echo "Total snps in {output.eva_filtered}: $(bcftools view --threads {threads} {output.eva_filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule ensemble_intersection:
    priority: 100
    input:
        bcftools = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz",
        bcftools_csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz.csi",
        freebayes = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.eva.vcf.gz",
        freebayes_csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.eva.vcf.gz.csi",
        haplotypecaller = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.eva.vcf.gz",
        haplotypecaller_csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.eva.vcf.gz.csi",
    output:
        bcftools_common = "results/05_ensemble/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.intersect.vcf.gz",
        freebayes_common = "results/05_ensemble/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.eva.intersect.vcf.gz",
        haplotypecaller_common = "results/05_ensemble/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.eva.intersect.vcf.gz",
    threads:6
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

        bcftools isec -O z8 -p results/05_ensemble -n~111 --threads {threads} {input.bcftools} {input.freebayes} {input.haplotypecaller};

        mv results/05_ensemble/0000.vcf.gz {output.bcftools_common};
        echo "Total snps in {output.bcftools_common}: $(bcftools view --threads {threads} {output.bcftools_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        mv results/05_ensemble/0001.vcf.gz {output.freebayes_common};
        echo "Total snps in {output.freebayes_common}: $(bcftools view --threads {threads} {output.freebayes_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        mv results/05_ensemble/0002.vcf.gz {output.haplotypecaller_common};
        echo "Total snps in {output.haplotypecaller_common}: $(bcftools view --threads {threads} {output.haplotypecaller_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        bcftools index --threads {threads} {output.bcftools_common};

        bcftools index --threads {threads} {output.freebayes_common};

        bcftools index --threads {threads} {output.haplotypecaller_common};

        rm results/05_ensemble/0000.vcf.gz.tbi
        rm results/05_ensemble/0001.vcf.gz.tbi
        rm results/05_ensemble/0002.vcf.gz.tbi

        """


rule view_bcftools_regions: #TODO Updated files
    priority:100
    input:
        filtered = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz",
        regions = "resources/snp_targets_100bp.txt"
    output:
        filtered_snps = "results/03_filtered/merged.filteredsnvs.QUAL60.bcftools.vcf.gz.pigmentSNPs.vcf"
    threads:2
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
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
    threads:2
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
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
    threads:2
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        " bcftools view "
        "-R {input.regions} "
        "--threads {threads} "
        "{input.filtered} "
        "> {output.filtered_snps} "


