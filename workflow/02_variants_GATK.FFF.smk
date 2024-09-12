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
CHROM = ('NC_056054.1', 'NC_056055.1', 'NC_056056.1', 'NC_056057.1', 'NC_056058.1', 'NC_056059.1', 'NC_056060.1', 'NC_056061.1', 'NC_056062.1', 'NC_056063.1', 'NC_056064.1', 'NC_056065.1', 'NC_056066.1', 'NC_056067.1', 'NC_056068.1', 'NC_056069.1', 'NC_056070.1', 'NC_056071.1', 'NC_056072.1', 'NC_056073.1', 'NC_056074.1', 'NC_056075.1', 'NC_056076.1', 'NC_056077.1', 'NC_056078.1', 'NC_056079.1', 'NC_056080.1')


rule all:
    input:
        "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi"



rule gatk_HaplotypeCaller_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna",
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
        partition = "compute",
        DTMP = "tmp",
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
    priority:1000
    input:
        vcfgz = expand("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz", chromosome = CHROM),
        csi = expand("results/02_snvs/{samples}.rawsnvs.{chromosome}.haplotypeCaller.vcf.gz", chromosome = CHROM),
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
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools concat --threads {threads} {input.vcfgz} | bcftools sort -O z8 -o {output.merged} - ;

        bcftools index --threads {threads} {output.merged} -o {output.csi}

        """


rule bcftools_norm_samples:
    priority: 100
    input:
        unnormal = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi",
    output:
        norm = temp("results/02_snvs/{samples}.rawsnvs.norm.haplotypeCaller.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.rawsnvs.norm.haplotypeCaller.vcf.gz.csi"),
    threads:6
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools norm --threads {threads} -O z8 -m- -f resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -o {output.norm} {input.unnormal};
    
        bcftools index --threads {threads} {output.norm};
    
        """


rule filter_DP:
    priority:100
    input:
        norm = "results/02_snvs/{samples}.rawsnvs.norm.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.norm.haplotypeCaller.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.haplotypeCaller.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.haplotypeCaller.vcf.gz.csi"),
    threads: 8
    conda:
        "bcftools"
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


rule filter_QUAL60: 
    priority:100
    input:
        dpfiltered = "results/03_filtered/{samples}.rawsnvs.norm.DPFilt.haplotypeCaller.vcf.gz",
        csi = "results/03_filtered/{samples}.rawsnvs.norm.DPFilt.haplotypeCaller.vcf.gz.csi",
    output:
        filtered = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz"),
        csi = temp("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools"
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


rule merge_animals_vcf:
    priority:100
    input:
        vcfgz = expand("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz", samples = SAMPLES),
        csi = expand("results/03_filtered/{samples}.rawsnvs.norm.DPFilt.QUAL60.bcftools.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi",
    benchmark:
        "benchmarks/merge_animals_vcf.tsv"
    threads: 16
    conda:
        "bcftools"
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


rule view_haplotype_chrom:
    priority:100
    input:
        merged_vcf = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/merged.FFF.rawsnvs.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi",
    output:
        filtered_vcf = "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz",
        filtered_vcf_csi = "results/02_snvs/merged.FFF.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi"
    params:
        chromosomes = "NC_056054.1,NC_056055.1,NC_056056.1,NC_056057.1,NC_056058.1,NC_056059.1,NC_056060.1,NC_056061.1,NC_056062.1,NC_056063.1,NC_056064.1,NC_056065.1,NC_056066.1,NC_056067.1,NC_056068.1,NC_056069.1,NC_056070.1,NC_056071.1,NC_056072.1,NC_056073.1,NC_056074.1,NC_056075.1,NC_056076.1,NC_056077.1,NC_056078.1,NC_056079.1,NC_056080.1"

    benchmark:
        "benchmarks/view_haplotypeCaller_chrom.tsv"
    threads: 8
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.merged_vcf} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes};

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi}; 

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """
