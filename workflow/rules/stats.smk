#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 stats.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile shave2.smk --cores X --use-conda
# Latest modification:  2023.02.02
# Done:                 Divided snakefile rules
###############################################################################
# CONFIGURATION FILES #
configfile: "config/config.yaml"
cluster_config: "slurm/config.yaml"

import os, sys
from snakemake.utils import min_version
min_version("5.18.0")

###############################################################################
# WILDCARDS #
SAMPLE, = glob_wildcards("resources/reads/{sample}_R1.fastq.gz")

###############################################################################
# RESOURCES #
TMPDIR = config["resources"]["tmpdir"] # Temporary directory

###############################################################################
# PARAMETERS #
reference_file = config["ref"]["reference"]

###############################################################################
rule samtools_stats:
    message:
        "SamTools stats: Collects statistics from BAM files"
    threads: 1
    resources: 
        partition='fast',
        mem_mb=4000,
        runtime=120,
    input:
        bam = rules.fixmateinformation.output.fixed,
        ref = reference_file,
    output:
        stats = "results/00_Quality_Control/{sample}_md_fixed_stats.txt",
    log:
        "results/11_Reports/samtools/{sample}_md_fixed_stats.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools stats --threads {threads} -r {input.ref} {input.bam} 1> {output.stats} 2> {log}
        """

###############################################################################
rule validate_sam:
    message:
        "Picard ValidateSamFile: Basic check for bam file validity, as interpreted by the Broad Institute."
    threads: 1
    resources: 
        partition='fast',
        mem_mb=4000, 
        runtime=600,
    input:
        bam = rules.fixmateinformation.output.fixed
    output:
        check = "results/00_Quality_Control/validatesamfile/{sample}_md_fixed_ValidateSam.txt",
    log:
        "results/11_Reports/validatesamfiles/{sample}_md_fixed_validate_bam.log",
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
        picard ValidateSamFile -I {input.bam} -O {output.check} -M SUMMARY > {log} 2>&1 || true
        """

###############################################################################
rule samtools_idxstats:
    message:
        "samtools idxstats: reports alignment summary statistics"
    threads: 1
    resources: 
        partition='fast',
        mem_mb=4000, 
        runtime=600,
    params:
        extra="",  # optional params string
    input:
        bam = rules.fixmateinformation.output.fixed,
        idx = rules.index_fixed_bam.output.index,
    output:
        idxstats = "results/00_Quality_Control/{sample}_md_fixed.idxstats.txt",
    log:
        "results/11_Reports/samtools/idxstats/{sample}_md_fixed_idxstats.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools idxstats {input.bam} > {output} &> {log}
        """

###############################################################################
rule samtools_flagstat:
    message:
        "samtools flagstats"
    threads: 1
    resources: 
        partition='fast',
        mem_mb=4000, 
        runtime=600,
    params:
        extra="",  # optional params string
    input:
        bam = rules.fixmateinformation.output.fixed,
    output:
        flagstat = "results/00_Quality_Control/{sample}_md_fixed_bam.flagstat.txt",
    log:
        "results/11_Reports/samtools/flagstat/{sample}_md_fixed_bam.log",

    shell:
        config["MODULES"]["SAMTOOLS"]+"""
        samtools flagstat {input.bam} > {output.flagstat} &> {log}
        """

###############################################################################
rule multiqc:
    message:
        "MultiQC"
    threads: 1
    resources: 
        partition='fast',
        mem_mb=8000, 
        runtime=1200,
    input:
        rules.samtools_flagstat.output.flagstat,
        rules.samtools_idxstats.output.idxstats,
        rules.samtools_stats.output.stats,
        rules.mark_duplicates.output.metrics,
        "results/00_Quality_Control/fastqc/",
        "results/00_Quality_Control/fastq-screen/",
    output:
        "multiqc_report.html"

    log:
        "results/11_Reports/multiqc/multiqc.log"
    shell:
        config["MODULES"]["MULTIQC"]+"""
        multiqc {input} -o results/00_Quality_Control/MULTIQC/ -n {output} > {log} 2>&1 
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

###############################################################################
rule qualimap:
    message:
        "Qualimap"
    threads: 1
    resources: 
        partition='fast',
        mem_mb=8000,
        runtime=1200,
    params:
        outdir = "results/00_Quality_Control/qualimap/{sample}/"
    input:
        bam = rules.fixmateinformation.output.fixed,
    output:
        protected("results/00_Quality_Control/qualimap/{sample}/qualimapReport.html"),
        protected("results/00_Quality_Control/qualimap/{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}/genome_results.txt"),
        protected("results/00_Quality_Control/qualimap/{sample}/raw_data_qualimapReport/coverage_histogram.txt")

    log:
        stderr="results/11_Reports/qualimap/logs/{sample}_qualimap.stderr",
        stdout="results/11_Reports/qualimap/logs/{sample}_qualimap.stdout",
    shell:
        config["MODULES"]["QUALIMAP"]+"""
            unset DISPLAY && qualimap bamqc -bam {input.bam} -nt {threads} --java-mem-size=8G -outdir {params.outdir} 
            exitcode=$?
            if [ $exitcode -eq 1 ]
            then
                exit 1
            else
                exit 0
            fi
        """
