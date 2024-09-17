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


SAMPLES = ('1941', '1927', '1936', '1938')


rule all:
    input:
        "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.vcf.gz",
        "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.vcf.gz",


### bcftools 

rule bcftools_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna",
    output:
        vcf = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz", # removed temp
        csi = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz.csi", # removed temp
    log:
        "logs/bcftools_vcf.{samples}.log"
    benchmark:
        "benchmarks/bcftools_vcf.{samples}.tsv"
    threads: 24
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 720 + ((attempt - 1) * 720),
        partition = "compute",
        DTMP = "tmp",
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
        "| bcftools view --threads {threads} -O z8 -o {output.vcf}; "
        "bcftools index --threads {threads} {output.vcf} -o  {output.csi} "


rule view_bcftools_chrom: #TODO
    priority:100
    input:
        merged = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz",
        csi = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi"
    output:
        filtered_vcf = "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.vcf.gz",
        filtered_vcf_csi = "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi"
    params:
        chromosomes = "NC_056054.1,NC_056055.1,NC_056056.1,NC_056057.1,NC_056058.1,NC_056059.1,NC_056060.1,NC_056061.1,NC_056062.1,NC_056063.1,NC_056064.1,NC_056065.1,NC_056066.1,NC_056067.1,NC_056068.1,NC_056069.1,NC_056070.1,NC_056071.1,NC_056072.1,NC_056073.1,NC_056074.1,NC_056075.1,NC_056076.1,NC_056077.1,NC_056078.1,NC_056079.1,NC_056080.1"
    benchmark:
        "benchmarks/view_bcftools_chrom.tsv"
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.merged} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes} &&

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi} &&

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """


rule norm_samples_bcftools:
    priority: 100
    input:
        unnormal = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz.csi",
    output:
        norm = temp("results/02_snvs/{samples}.rawsnvs.norm.bcftools.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.norm.bcftools.vcf.gz.csi"),
    threads:6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools norm --threads {threads} -O z8 -m+ -s -D -f resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -o {output.norm} {input.unnormal};
    
        bcftools index --threads {threads} {output.norm};
    
        """


rule filter_DP_bcftools:
    priority:100
    input:
        norm = "results/02_snvs/{samples}.rawsnvs.norm.bcftools.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.norm.bcftools.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.bcftools.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.bcftools.vcf.gz.csi"),
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        # -e is 'exclude'

        bcftools view --threads {threads} -e 'INFO/DP<5 || INFO/DP>250' {input.norm} -O z8 -o {output.filtered} - ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule filter_QUAL60_bcftools: 
    priority:100
    input:
        dpfiltered = "results/03_filtered/{samples}.rawsnvs.norm.DPFilt.bcftools.vcf.gz",
        csi = "results/03_filtered/{samples}.rawsnvs.norm.DPFilt.bcftools.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''
        # -e is 'exclude'

        bcftools view -e 'QUAL<60' {input.dpfiltered} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        '''


rule merge_bcftools_vcf:
    priority:100
    input:
        vcfgz = expand("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz", samples = SAMPLES),
        csi = expand("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz",
        csi = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi"
    benchmark:
        "benchmarks/merge_bcftools_vcf.tsv"
    threads: 16
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
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
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna",
    output:
        vcfgz = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz", # removed temp
        csi = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz.csi", # removed temp
    log:
        "logs/freebayes_vcf.{samples}.log"
    benchmark:
        "benchmarks/freebayes_vcf.{samples}.tsv"
    threads: 2
    conda:
        "freebayes-1.3.8"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 12),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "freebayes "
        "--standard-filters " # Use stringent input base and mapping quality filters. Equivalent to -m 30 -q 20 -R 0 -S 0
        "--limit-coverage 250 " # Match the other variant callers
        #"--pooled-continuous " # Output all alleles which pass input filters, regardles of genotyping outcome or model.
        #"--trim-complex-tail " # Trim complex tails.
        #"-F 0.01 " # minimum fraction of observations supporting alternate allele within one individual [0.05]
        "-f {input.referenceGenome} {input.bam} "
        "| bcftools view --threads {threads} -O z8 -o {output.vcfgz}; "
        "bcftools index --threads {threads} {output.vcfgz} -o  {output.csi} "


rule view_freebayes_chrom: #TODO
    priority:100
    input:
        merged = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz",
        csi = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz.csi"
    output:
        filtered_vcf = "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.vcf.gz",
        filtered_vcf_csi = "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.freebayes.vcf.gz.csi"
    params:
        chromosomes = "NC_056054.1,NC_056055.1,NC_056056.1,NC_056057.1,NC_056058.1,NC_056059.1,NC_056060.1,NC_056061.1,NC_056062.1,NC_056063.1,NC_056064.1,NC_056065.1,NC_056066.1,NC_056067.1,NC_056068.1,NC_056069.1,NC_056070.1,NC_056071.1,NC_056072.1,NC_056073.1,NC_056074.1,NC_056075.1,NC_056076.1,NC_056077.1,NC_056078.1,NC_056079.1,NC_056080.1"
    benchmark:
        "benchmarks/view_freebayes_chrom.tsv"
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.merged} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes} &&

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi} &&

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """


rule norm_samples_freebayes:
    priority: 100
    input:
        unnormal = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf.gz.csi",
    output:
        norm = temp("results/02_snvs/{samples}.rawsnvs.norm.freebayes.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.norm.freebayes.vcf.gz.csi"),
    threads:6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools norm --threads {threads} -O z8 -m+ -s -D -f resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -o {output.norm} {input.unnormal};
    
        echo "Total snps in {output.norm}: $(bcftools view --threads {threads} {output.norm} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        bcftools index --threads {threads} {output.norm};
    
        """


rule filter_DP_freebayes:
    priority:100
    input:
        norm = "results/02_snvs/{samples}.rawsnvs.norm.freebayes.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.norm.freebayes.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.freebayes.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.freebayes.vcf.gz.csi"),
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        # -e is 'exclude'

        bcftools view --threads {threads} -e 'INFO/DP<5 || INFO/DP>250' {input.norm} -O z8 -o {output.filtered} - ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule filter_QUAL60_freebayes: 
    priority:100
    input:
        dpfiltered = "results/03_filtered/{samples}.rawsnvs.norm.DPFilt.freebayes.vcf.gz",
        csi = "results/03_filtered/{samples}.rawsnvs.norm.DPFilt.freebayes.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''
        # -e is 'exclude'

        bcftools view -e 'QUAL<60' {input.dpfiltered} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        '''


rule merge_freebayes_vcf:
    priority:100
    input:
        vcfgz = expand("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz", samples = SAMPLES),
        csi = expand("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz",
        csi = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.freebayes.vcf.gz.csi"
    benchmark:
        "benchmarks/merge_freebayes_vcf.tsv"
    threads: 16
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools merge --threads {threads} {input.vcfgz} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 


        """

