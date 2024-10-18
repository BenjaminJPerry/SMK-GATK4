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
        "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz",


rule get_eva_snvs:
    output:
        vcf = "resources/eva/9913_GCA_002263795.2_current_ids.vcf.gz",
        csi = "resources/eva/9913_GCA_002263795.2_current_ids.vcf.gz.csi",
    params:
        vcf = "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/bos_taurus/ARSUCD1.2/9913_GCA_002263795.2_current_ids.vcf.gz",
        index = "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/bos_taurus/ARSUCD1.2/9913_GCA_002263795.2_current_ids.vcf.gz.csi",
    threads: 2
    shell:
        """
        
        wget -t 5 -O {output.vcf} {params.vcf} ;
        wget -t 5 -O {output.csi} {params.index} ; 
        
        """


rule rename_eva_snvs:
    input:
        vcf = "resources/eva/9913_GCA_002263795.2_current_ids.vcf.gz",
        csi = "resources/eva/9913_GCA_002263795.2_current_ids.vcf.gz.csi",
        sed_file = "eva_rename.sed",
    output:
        vcf = "resources/eva/9913_GCA_002263795.2_current_ids.sed.vcf.gz",
        csi = "resources/eva/9913_GCA_002263795.2_current_ids.sed.vcf.gz.csi",
    threads: 8
    conda:
        "bcftools"
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
        eva = "resources/eva/9913_GCA_002263795.2_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9913_GCA_002263795.2_current_ids.sed.vcf.gz.csi"
    output:
        eva_filtered = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz",
        csi = "results/04_merged/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.eva.vcf.gz.csi",
    threads:8
    conda:
        "bcftools"
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


rule isec_freebayes_LIC565:
    priority:100
    input:
        filtered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz.csi",
        snps_LIC = "resources/LIC_565.ch.frmt.sorted.vcf.gz",
        snps_LIC_csi = "resources/LIC_565.ch.frmt.sorted.vcf.gz.csi",
    output:
        LICFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.vcf.gz.csi"
    threads:8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/03_filtered/isec_freebayes_LIC565 {input.snps_LIC} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_freebayes_LIC565/0001.vcf -o {output.LICFiltered} ;

        bcftools index --threads {threads} {output.LICFiltered} -o {output.csi} ; 

        rm -r results/03_filtered/isec_freebayes_LIC565 ;

        echo "Total snps in {output.LICFiltered}: $(bcftools view --threads {threads} {output.LICFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule isec_haplotypeCaller_LIC565:
    priority:100
    input:
        filtered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi",
        snps_LIC = "resources/LIC_565.ch.frmt.sorted.vcf.gz",
        snps_LIC_csi = "resources/LIC_565.ch.frmt.sorted.vcf.gz.csi",
    output:
        LICFiltered = temp("results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.vcf.gz"),
        csi = temp("results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.vcf.gz.csi")
    threads:8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/03_filtered/isec_haplotypeCaller_LIC565 {input.snps_LIC} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_haplotypeCaller_LIC565/0001.vcf -o {output.LICFiltered} ;

        bcftools index --threads {threads} {output.LICFiltered} -o {output.csi} ; 

        rm -r results/03_filtered/isec_haplotypeCaller_LIC565 ;

        echo "Total snps in {output.LICFiltered}: $(bcftools view --threads {threads} {output.LICFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule isec_bcftools_TBulls:
    priority:100
    input:
        LICFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.vcf.gz.csi",
        snps_TBulls = "resources/ARS_1000bullsgenome.frmt.sorted.vcf.gz",
        snps_TBulls_csi = "resources/ARS_1000bullsgenome.frmt.sorted.vcf.gz.csi",
    output:
        TBullsFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz.csi"
    threads:8
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
        bcftools isec -O v -p results/03_filtered/isec_bcftools_TBulls {input.snps_TBulls} {input.LICFiltered} ;

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_bcftools_TBulls/0001.vcf -o {output.TBullsFiltered} ;

        bcftools index --threads {threads} {output.TBullsFiltered} -o {output.csi} ; 

        rm -r results/03_filtered/isec_bcftools_TBulls ;

        echo "Total snps in {output.TBullsFiltered}: $(bcftools view --threads {threads} {output.TBullsFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule isec_freebayes_TBulls:
    priority:100
    input:
        LICFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.vcf.gz.csi",
        snps_TBulls = "resources/ARS_1000bullsgenome.frmt.sorted.vcf.gz",
        snps_TBulls_csi = "resources/ARS_1000bullsgenome.frmt.sorted.vcf.gz.csi",
    output:
        TBullsFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz.csi"
    threads:8
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
        bcftools isec -O v -p results/03_filtered/isec_freebayes_TBulls {input.snps_TBulls} {input.LICFiltered} ;

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_freebayes_TBulls/0001.vcf -o {output.TBullsFiltered} ;

        bcftools index --threads {threads} {output.TBullsFiltered} -o {output.csi} ; 

        rm -r results/03_filtered/isec_freebayes_TBulls ;

        echo "Total snps in {output.TBullsFiltered}: $(bcftools view --threads {threads} {output.TBullsFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule isec_haplotypeCaller_TBulls:
    priority:100
    input:
        LICFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.vcf.gz.csi",
        snps_TBulls = "resources/ARS_1000bullsgenome.frmt.sorted.vcf.gz",
        snps_TBulls_csi = "resources/ARS_1000bullsgenome.frmt.sorted.vcf.gz.csi",
    output:
        TBullsFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz.csi"
    threads:8
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
        bcftools isec -O v -p results/03_filtered/isec_haplotypeCaller_TBulls {input.snps_TBulls} {input.LICFiltered} ;

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_haplotypeCaller_TBulls/0001.vcf -o {output.TBullsFiltered} ;

        bcftools index --threads {threads} {output.TBullsFiltered} -o {output.csi} ; 

        rm -r results/03_filtered/isec_haplotypeCaller_TBulls ;

        echo "Total snps in {output.TBullsFiltered}: $(bcftools view --threads {threads} {output.TBullsFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule ensemble_intersection:
    priority: 100
    input:
        bcftools = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.vcf.gz",
        freebayes = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.vcf.gz",
        haplotypeCaller = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.norm.vcf.gz",
    output:
        bcftools_common = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.intersect.vcf.gz",
        freebayes_common = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.intersect.vcf.gz",
        haplotypeCaller_common = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.norm.intersect.vcf.gz",
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

        bcftools isec -O z8 -p results/05_ensemble -n~111 --threads {threads} {input.bcftools} {input.freebayes} {input.haplotypeCaller};

        mv results/05_ensemble/0000.vcf.gz {output.bcftools_common};
        echo "Total snps in {output.bcftools_common}: $(bcftools view --threads {threads} {output.bcftools_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        mv results/05_ensemble/0001.vcf.gz {output.freebayes_common};
        echo "Total snps in {output.freebayes_common}: $(bcftools view --threads {threads} {output.freebayes_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        mv results/05_ensemble/0002.vcf.gz {output.haplotypeCaller_common};
        echo "Total snps in {output.haplotypeCaller_common}: $(bcftools view --threads {threads} {output.haplotypeCaller_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        bcftools index --threads {threads} {output.bcftools_common};

        bcftools index --threads {threads} {output.freebayes_common};

        bcftools index --threads {threads} {output.haplotypeCaller_common};

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


