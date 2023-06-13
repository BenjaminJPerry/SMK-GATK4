# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"


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


SAMPLES, = glob_wildcards("fastq/AHJNKHDSX5{samples}_R1.fastq.gz")


rule all:
    input:
        "",

rule gatk_HaplotypeCaller_first_pass:
    input:
        bams = "results/mapped/{sample}_recalibrated.bam",
        refgenome = expand("{refgenome}", refgenome = config["REFGENOME"]),
    output:
        protected("results/called/{sample}_raw_snps_indels.vcf")
    params:
        tdir = config["TEMPDIR"],
        padding = get_wes_padding_command,
        intervals = get_wes_intervals_command
    log:
        "logs/gatk_HaplotypeCaller_single/{sample}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_single/{sample}.tsv"
    singularity:
        "docker://broadinstitute/gatk:4.2.6.1"
    threads: 2
    resources:
        mem_gb = 8,
        # partition = config["PARTITION"]["CPU"]
    message:
        "Calling germline SNPs and indels via local re-assembly of haplotypes for {input.bams}"
    shell:
        "gatk HaplotypeCaller "
        "--java-options '-Xmx{resources.mem_gb}G' "
        "-I {input.bams} "
        "-R {input.refgenome} "
        # "-D {input.dbsnp} "
        "-O {output} "
        "--tmp-dir {params.tdir} {params.padding} {params.intervals} "
        "&> {log}"


