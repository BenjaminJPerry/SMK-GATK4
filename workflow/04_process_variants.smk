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
#         "results/02_snvs/merged.chrom.bcftools.vcf.gz",
#         "results/02_snvs/merged.chrom.freebayes.vcf.gz",
#         "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz"

rule all:
    input:
        expand("results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz", samples = SAMPLES),
        expand("results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz", samples = SAMPLES),
        expand("results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz", samples = SAMPLES),
        "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz",
        "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz",
        "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.unique.vcf.gz",
        "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.unique.vcf.gz",
        "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypCaller.LIC565.TBulls.norm.unique.vcf.gz",
        "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.intersect.vcf.gz",
        "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.intersect.vcf.gz",
        "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypCaller.LIC565.TBulls.norm.intersect.vcf.gz",

        # #"results/03_filtered/merged.chrom.freebayes.QUAL60.vcf.gz.pigmentSNPs.vcf",
        #"results/03_filtered/merged.chrom.bcftools.QUAL60.vcf.gz.pigmentSNPs.vcf",
        #"results/03_filtered/merged.chrom.haplotypeCaller.QUAL60.vcf.gz.pigmentSNPs.vcf",


rule bcftools_private_snps:
    priority:100
    input:
        merged = "results/02_snvs/merged.chrom.bcftools.vcf.gz",
        csi = "results/02_snvs/merged.chrom.bcftools.vcf.gz.csi",
    output:
        private = temp("results/03_filtered/{samples}.chrom.private.bcftools.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.bcftools.vcf.gz.csi"),
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
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.merged} -o {output.private} ;

        bcftools index --threads {threads} {output.private} -o {output.csi}

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a bcftools.animals.private.snps.counts.summary.txt &&

        exit 0;

        """


rule freebayes_private_snps:
    priority:100
    input:
        merged = "results/02_snvs/merged.chrom.freebayes.vcf.gz",
        csi = "results/02_snvs/merged.chrom.freebayes.vcf.gz.csi",
    output:
        private = temp("results/03_filtered/{samples}.chrom.private.freebayes.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.freebayes.vcf.gz.csi"),
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
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.merged} -o {output.private} ;

        bcftools index --threads {threads} {output.private} -o {output.csi}

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a freebayes.animals.private.snps.counts.summary.txt &&

        exit 0;

        """


rule haplotypeCaller_private_snps:
    priority:100
    input:
        merged = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/merged.chrom.haplotypeCaller.vcf.gz.csi",
    output:
        private = temp("results/03_filtered/{samples}.chrom.private.haplotypeCaller.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.haplotypeCaller.vcf.gz.csi"),
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
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.merged} -o {output.private} ;

        bcftools index --threads {threads} {output.private} -o {output.csi}

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a haplotypeCaller.animals.private.snps.counts.summary.txt &&

        exit 0;

        """


rule bcftools_DPFilter_private:
    priority:100
    input:
        private = "results/03_filtered/{samples}.chrom.private.bcftools.vcf.gz",
        csi = "results/03_filtered/{samples}.chrom.private.bcftools.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.chrom.private.DPFilt.bcftools.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.DPFilt.bcftools.vcf.gz.csi"),
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

        bcftools norm --threads {threads} -m- -f /nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa {input.private} | bcftools view --threads {threads} -O z8 -e 'INFO/DP<10 || INFO/DP>2500' -o {output.filtered} - ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi}

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a bcftools.animals.private.snps.counts.summary.txt &&

        """


rule freebayes_DPFilter_private:
    priority:100
    input:
        private = "results/03_filtered/{samples}.chrom.private.freebayes.vcf.gz",
        csi = "results/03_filtered/{samples}.chrom.private.freebayes.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.chrom.private.DPFilt.freebayes.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.DPFilt.freebayes.vcf.gz.csi"),
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

        bcftools norm --threads {threads} -m- -f /nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa {input.private} | bcftools view --threads {threads} -O z8 -e 'INFO/DP<10 || INFO/DP>2500' -o {output.filtered} - ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi}
       
        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a freebayes.animals.private.snps.counts.summary.txt &&

        """


rule haplotypeCaller_DPFilter_private:
    priority:100
    input:
        private = "results/03_filtered/{samples}.chrom.private.haplotypeCaller.vcf.gz",
        csi = "results/03_filtered/{samples}.chrom.private.haplotypeCaller.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.chrom.private.DPFilt.haplotypeCaller.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.DPFilt.haplotypeCaller.vcf.gz.csi"),
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

        bcftools norm --threads {threads} -m- -f /nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa {input.private} | bcftools view --threads {threads} -O z8 -e 'INFO/DP<10 || INFO/DP>2500' -o {output.filtered} - ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi}

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a haplotypeCaller.animals.private.snps.counts.summary.txt &&

        """


rule bcftools_filter_private_QUAL60: 
    priority:100
    input:
        DPFilt = "results/03_filtered/{samples}.chrom.private.DPFilt.bcftools.vcf.gz",
        index = "results/03_filtered/{samples}.chrom.private.DPFilt.bcftools.vcf.gz.csi"
    output:
        filtered = temp("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz.csi"),
    threads:8
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

        bcftools view -e 'QUAL<60' {input.DPFilt} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi} 

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a bcftools.animals.private.snps.counts.summary.txt &&

        '''


rule freebayes_filter_private_QUAL60:
    priority:100
    input:
        DPFilt = "results/03_filtered/{samples}.chrom.private.DPFilt.freebayes.vcf.gz",
        index = "results/03_filtered/{samples}.chrom.private.DPFilt.freebayes.vcf.gz.csi"
    output:
        filtered = temp("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz.csi"),
    threads:8
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
        
        bcftools view -e 'QUAL<60' {input.DPFilt} -O z8 -o {output.filtered};
        
        bcftools index --threads {threads} {output.filtered} -o {output.csi} 

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a freebayes.animals.private.snps.counts.summary.txt &&

        '''


rule haplotypeCaller_filter_private_QUAL60:
    priority:100
    input:
        DPFilt = "results/03_filtered/{samples}.chrom.private.DPFilt.haplotypeCaller.vcf.gz",
        index = "results/03_filtered/{samples}.chrom.private.DPFilt.haplotypeCaller.vcf.gz.csi"
    output:
        filtered = temp("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi"),
    threads:8
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
        
        bcftools view -e 'QUAL<60' {input.DPFilt} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi} 

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a haplotypeCaller.animals.private.snps.counts.summary.txt &&

        '''


rule bcftools_gather_private:
    priority:100
    input:
        filtered = expand("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz", samples=SAMPLES),
        csi = expand("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz.csi", samples=SAMPLES),
    output:
        merged = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz.csi",
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
        bcftools merge --threads {threads} {input.filtered} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;
        echo "Total snps in {output.merged} at QUAL>=60: $(bcftools view --threads {threads} -i 'QUAL>=60' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=50: $(bcftools view --threads {threads} -i 'QUAL>=50' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=40: $(bcftools view --threads {threads} -i 'QUAL>=40' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=30: $(bcftools view --threads {threads} -i 'QUAL>=30' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=20: $(bcftools view --threads {threads} -i 'QUAL>=20' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt &&

        """


rule freebayes_gather_private:
    priority:100
    input:
        filtered = expand("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz", samples=SAMPLES),
        csi = expand("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz.csi", samples=SAMPLES),
    output:
        merged = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.vcf.gz.csi",
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
        bcftools merge --threads {threads} {input.filtered} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;
        echo "Total snps in {output.merged} at QUAL>=60: $(bcftools view --threads {threads} -i 'QUAL>=60' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=50: $(bcftools view --threads {threads} -i 'QUAL>=50' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=40: $(bcftools view --threads {threads} -i 'QUAL>=40' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=30: $(bcftools view --threads {threads} -i 'QUAL>=30' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=20: $(bcftools view --threads {threads} -i 'QUAL>=20' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt &&

        """



rule haplotypeCaller_gather_private:
    priority:100
    input:
        filtered = expand("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz", samples=SAMPLES),
        csi = expand("results/03_filtered/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi", samples=SAMPLES),
    output:
        merged = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi",
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
        bcftools merge --threads {threads} {input.filtered} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;
        echo "Total snps in {output.merged} at QUAL>=60: $(bcftools view --threads {threads} -i 'QUAL>=60' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=50: $(bcftools view --threads {threads} -i 'QUAL>=50' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=40: $(bcftools view --threads {threads} -i 'QUAL>=40' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=30: $(bcftools view --threads {threads} -i 'QUAL>=30' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        echo "Total snps in {output.merged} at QUAL>=20: $(bcftools view --threads {threads} -i 'QUAL>=20' {output.merged} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt &&


        """


rule isec_bcftools_LIC565:
    priority:100
    input:
        filtered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.vcf.gz.csi",
        snps_LIC = "resources/LIC_565.ch.frmt.sorted.vcf.gz",
        snps_LIC_csi = "resources/LIC_565.ch.frmt.sorted.vcf.gz.csi",
    output:
        LICFiltered = temp("results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.vcf.gz"),
        csi = temp("results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.vcf.gz.csi"),
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
        bcftools isec -O v -p results/03_filtered/isec_bcftools_LIC565 {input.snps_LIC} {input.filtered} &&

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_bcftools_LIC565/0001.vcf -o {output.LICFiltered} &&

        bcftools index --threads {threads} {output.LICFiltered} -o {output.csi} && 

        rm -r results/03_filtered/isec_bcftools_LIC565 &&

        echo "Total snps in {output.LICFiltered}: $(bcftools view --threads {threads} {output.LICFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

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
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/03_filtered/isec_freebayes_LIC565 {input.snps_LIC} {input.filtered} &&

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_freebayes_LIC565/0001.vcf -o {output.LICFiltered} &&

        bcftools index --threads {threads} {output.LICFiltered} -o {output.csi} && 

        rm -r results/03_filtered/isec_freebayes_LIC565 &&

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
        partition = "large,milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/03_filtered/isec_haplotypeCaller_LIC565 {input.snps_LIC} {input.filtered} &&

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_haplotypeCaller_LIC565/0001.vcf -o {output.LICFiltered} &&

        bcftools index --threads {threads} {output.LICFiltered} -o {output.csi} && 

        rm -r results/03_filtered/isec_haplotypeCaller_LIC565 &&

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
        bcftools isec -O v -p results/03_filtered/isec_bcftools_TBulls {input.snps_TBulls} {input.LICFiltered} &&

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_bcftools_TBulls/0001.vcf -o {output.TBullsFiltered} &&

        bcftools index --threads {threads} {output.TBullsFiltered} -o {output.csi} && 

        rm -r results/03_filtered/isec_bcftools_TBulls &&

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
        bcftools isec -O v -p results/03_filtered/isec_freebayes_TBulls {input.snps_TBulls} {input.LICFiltered} &&

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_freebayes_TBulls/0001.vcf -o {output.TBullsFiltered} &&

        bcftools index --threads {threads} {output.TBullsFiltered} -o {output.csi} && 

        rm -r results/03_filtered/isec_freebayes_TBulls &&

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
        bcftools isec -O v -p results/03_filtered/isec_haplotypeCaller_TBulls {input.snps_TBulls} {input.LICFiltered} &&

        bcftools view --threads {threads} -O z8 results/03_filtered/isec_haplotypeCaller_TBulls/0001.vcf -o {output.TBullsFiltered} &&

        bcftools index --threads {threads} {output.TBullsFiltered} -o {output.csi} && 

        rm -r results/03_filtered/isec_haplotypeCaller_TBulls &&

        echo "Total snps in {output.TBullsFiltered}: $(bcftools view --threads {threads} {output.TBullsFiltered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


rule bcftools_final_private_snps:
    priority:100
    input:
        TBullsFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz.csi",
    output:
        private = "results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz",
        csi = "results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz.csi",
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
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.TBullsFiltered} -o {output.private} ;

        bcftools index --threads {threads} {output.private} -o {output.csi}

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a bcftools.animals.private.snps.counts.summary.txt &&

        exit 0;

        """


rule freebayes_final_private_snps:
    priority:100
    input:
        TBullsFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz.csi",
    output:
        private = "results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz",
        csi = "results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz.csi",
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
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.TBullsFiltered} -o {output.private} ;

        bcftools index --threads {threads} {output.private} -o {output.csi}

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a freebayes.animals.private.snps.counts.summary.txt &&

        exit 0;

        """


rule haplotypeCaller_final_private_snps:
    priority:100
    input:
        TBullsFiltered = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz",
        csi = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz.csi",
    output:
        private = "results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz",
        csi = "results/04_animals/{samples}.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz.csi",
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
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.TBullsFiltered} -o {output.private} ;

        bcftools index --threads {threads} {output.private} -o {output.csi}

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a haplotypeCaller.animals.private.snps.counts.summary.txt &&

        exit 0;

        """


rule bcftools_norm_merged:
    priority: 100
    input:
        unnormal = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.vcf.gz",
    output:
        normal = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.vcf.gz",
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

    bcftools norm --threads {threads} -O z8 -m- -f /nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa -o {output.normal} {input.unnormal};
    
    bcftools index --threads {threads} {output.normal};
    
    """


rule freebayes_norm_merged:
    priority: 100
    input:
        unnormal = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.vcf.gz",
    output:
        normal = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.vcf.gz",
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

    bcftools norm --threads {threads} -O z8 -m- -f /nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa -o {output.normal} {input.unnormal};
    
    bcftools index --threads {threads} {output.normal};
    
    """


rule haplotypeCaller_norm_merged:
    priority: 100
    input:
        unnormal = "results/03_filtered/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.vcf.gz",
    output:
        normal = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.norm.vcf.gz",
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

    bcftools norm --threads {threads} -O z8 -m- -f /nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa -o {output.normal} {input.unnormal};
    
    bcftools index --threads {threads} {output.normal};
    
    """


rule ensemble_intersection:
    priority: 100
    input:
        bcftools = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.vcf.gz",
        freebayes = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.vcf.gz",
        haplotypeCaller = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.norm.vcf.gz",
    output:=
        bcftools_uniq = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.bcftools.LIC565.TBulls.norm.unique.vcf.gz",
        freebayes_uniq = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.freebayes.LIC565.TBulls.norm.unique.vcf.gz",
        haplotypeCaller_uniq = "results/05_ensemble/merged.chrom.private.DPFilt.QUAL60.haplotypeCaller.LIC565.TBulls.norm.unique.vcf.gz",
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

        bcftools isec -O z8 -p results/05_ensemble -n=3 --threads {threads} {input.bcftools} {input.freebayes} {input.haplotypeCaller};

        mv 0001.vcf.gz {output.bcftools_uniq}
        echo "Total snps in {output.bcftools_uniq}: $(bcftools view --threads {threads} {output.bcftools_uniq} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt

        mv 0002.vcf.gz {output.freebayes_uniq}
        echo "Total snps in {output.freebayes_uniq}: $(bcftools view --threads {threads} {output.freebayes_uniq} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt

        mv 0003.vcf.gz {output.haplotypeCaller_uniq}
        echo "Total snps in {output.haplotypeCaller_uniq}: $(bcftools view --threads {threads} {output.haplotypeCaller_uniq} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt

        bcftools view -O z8 --threads {threads} -R results/05_ensemble/sites.txt {input.bcftools} -o {output.bcftools_common}
        bcftools index --threads {threads} {output.bcftools_common}
        echo "Total snps in {output.bcftools_common}: $(bcftools view --threads {threads} {output.bcftools_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt

        bcftools view -O z8 --threads {threads} -R results/05_ensemble/sites.txt {input.freebayes} -o {output.freebayes_common}
        bcftools index --threads {threads} {output.freebayes_common}
        echo "Total snps in {output.freebayes_common}: $(bcftools view --threads {threads} {output.freebayes_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt

        bcftools view -O z8 --threads {threads} -R results/05_ensemble/sites.txt {input.haplotypeCaller} -o {output.haplotypeCaller_common}
        bcftools index --threads {threads} {output.haplotypeCaller_common}
        echo "Total snps in {output.haplotypeCaller_common}: $(bcftools view --threads {threads} {output.haplotypeCaller_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt

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


