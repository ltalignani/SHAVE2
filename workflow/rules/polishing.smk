#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 polishing.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile polishing.smk --cores X --use-conda
# Latest modification:  2023.02.18
# Done:                 
###############################################################################
# CONFIGURATION FILES #
configfile: "config/config.yaml"
cluster_config: "slurm/config.yaml"

import os, sys
from snakemake.utils import min_version
min_version("6.12.0")

###############################################################################
# PARAMETERS #
reference_file = config["ref"]["reference"]

############################### R U L E S #####################################
rule SetNmMdAndUqTags:
    message:
        "Picard SetNmMdAndUqTags: this tool takes in a coordinate-sorted SAM or BAM and calculates the NM, MD, and UQ tags by comparing with the reference."
    resources: 
        partition='fast',
        mem_mb=4000,
        runtime=120,
    input:
        bam = rules.mark_duplicates.output.bam,
        ref = reference_file
    output:             
        fix = temp("results/04_Polishing/{sample}_sorted-mark-dup-fx.bam"),
    log:
        "results/11_Reports/SetNmMdAndUqTags/{sample}_sorted-mark-dup-fx.log",
    benchmark:
        "benchmarks/setnmmdanduqtags/{sample}.tsv",
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
            picard SetNmMdAndUqTags R={input.ref} I={input.bam} O={output.fix} > {log} 2>&1
        """

###############################################################################
rule fixmateinformation:
    message:
        "Picard FixMateInformation: This tool ensures that all mate-pair information is in sync between each read and its mate pair."
    resources: 
        partition='fast',
        mem_mb=4000,
        runtime=1200,
    input:
        bam = rules.SetNmMdAndUqTags.output.fix
    output:
        fixed = "results/04_Polishing/{sample}_md_fixed.bam",
    log:
        "results/11_Reports/fixmateinformation/{sample}_bwa_fixed.log",
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
            picard FixMateInformation -I {input.bam} -O {output.fixed} --ADD_MATE_CIGAR true &> {log}
        """

###############################################################################
rule index_fixed_bam:
    message:
        "SamTools indexing indel qualities BAM file {wildcards.sample} sample",
    resources: 
        partition='fast',
        cpus_per_task=4,
        mem_mb=4000,
        runtime=1200,
    input:
        fixmate = rules.fixmateinformation.output.fixed
    output:
        index   = "results/04_Polishing/{sample}_md_fixed.bai",
    log:
        "results/11_Reports/samtools/{sample}_fixed-mate-index.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
            samtools index -@ {resources.cpus_per_task} -b {input.fixmate} {output.index} &> {log}
        """

###############################################################################

