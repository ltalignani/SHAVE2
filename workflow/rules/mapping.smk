#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 mapping.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile workflow/rules/mapping.smk --cores X 
# Latest modification:  2023.07.13
# Done:                 Updated rules
###############################################################################
# CONFIGURATION FILES #
configfile: "config/config.yaml"
#cluster_config: "slurm/config.yaml"

import os, glob, sys, re
from snakemake.utils import min_version
min_version("6.12.0")

###############################################################################
# WILDCARDS #
# SAMPLES, READS, = glob_wildcards('raw/{sample}_R{read}.fastq')
SAMPLE, = glob_wildcards("raw/{sample}_R1.fastq.gz")
READS = ['1', '2']

###############################################################################
# RESOURCES #
TMPDIR = config["resources"]["tmpdir"] # Temporary directory

###############################################################################
# PARAMETERS #
CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset
REFPATH = config["ref"]["path"]                     # Path to genomes references
INDEXPATH = config["bwa"]["path"]                   # Genome reference index path
reference_file =  config["ref"]["reference"]
basename_reference = Path(reference_file).stem # Pathlib PurePath.stem: remove suffix
INDEX = config["ref"]["ref_name"]                   # Genome reference index files ['amb','ann','bwt','pac','sa']

###############################################################################
# FUNCTIONS AND COMMANDS #

############################## O N S T A R T ###################################

############################### R U L E S #####################################
rule fastqc_quality_control:
    message: "FastQC reads quality controling"
    resources: 
        partition='fast',
        mem_mb=4000,
        runtime=60,
    input:
        os.path.join("raw/", '{sample}_R{read}.fastq.gz')
    output:
        fastqc = os.path.join("results/00_Quality_Control/fastqc/", '{sample}_R{read}_fastqc.html'),
    log:
        "results/11_Reports/quality/{sample}_R{read}_fastqc.log"
    shell:
        config["MODULES"]["FASTQC"]+"""
            fastqc --quiet --threads {resources.cpus_per_task} --outdir {output.fastqc} {input} &> {log}
        """
###############################################################################
rule fastqscreen_contamination_checking:
    message: "Fastq-Screen reads contamination checking"
    resources: 
        partition='fast',
        mem_mb=4000,
        runtime=60,
    params:
        config = CONFIG,
        mapper = MAPPER,
        subset = SUBSET
    input:
        fastq = "raw/"
    output:
        fastqscreen = directory("results/00_Quality_Control/fastq-screen/"),
    log:
        "results/11_Reports/quality/fastq-screen.log"
    shell:
        config["MODULES"]["FASTQSCREEN"]+"\n"+config["MODULES"]["BWA"]+"""
            fastq_screen -q --threads {resources.cpus_per_task} --conf {params.config} --aligner {params.mapper} --subset {params.subset} --outdir {output.fastqscreen} {input.fastq}/*.fastq.gz &> {log}
        """

###############################################################################
rule trimmomatic:
    message: "Trimming reads for {wildcards.sample}"
    resources:
        partition='long',
        cpus_per_task=8,
        mem_mb=6000,
        runtime=300,
    input:
        r1="raw/{sample}_R1.fastq.gz",
        r2="raw/{sample}_R2.fastq.gz",
        adapters = config["trimmomatic"]["adapters"]["truseq2-pe"]
    output:
        forward_reads   = "results/01_Trimming/{sample}_trimmomatic_R1.fastq.gz",
        reverse_reads   = "results/01_Trimming/{sample}_trimmomatic_R2.fastq.gz",
        forwardUnpaired = "results/01_Trimming/{sample}_trimmomatic_unpaired_R1.fastq.gz",
        reverseUnpaired = "results/01_Trimming/{sample}_trimmomatic_unpaired_R2.fastq.gz"
    log:
        "results/11_Reports/trimmomatic/{sample}.log"
    params:
        settings = config["trimmomatic"]["settings"],
        phred    = config["trimmomatic"]["phred"]
    shell:
        config["MODULES"]["TRIMMOMATIC"]+"""
            trimmomatic PE -threads {resources.cpus_per_task} {params.phred} -trimlog {log} {input.r1} {input.r2} {output.forward_reads} {output.forwardUnpaired} {output.reverse_reads} {output.reverseUnpaired} {params.settings}
        """

###############################################################################
rule bwa_mapping:
    message:
        "BWA-MEM mapping sample reads against reference genome sequence"
    resources: 
        partition='long',
        cpus_per_task=16,
        mem_mb=32000,
        runtime=1200,
        slurm_extra="--mail-type=ALL --mail-user=loic.talignani@ird.fr"
    params: 
        #rg = f"@RG\\tID:{{sample}}\\tSM:{{sample}}\\tPL:Illumina", # Manage ReadGroup
        ref = reference_file,
        other_options_samtools_view = "-bh",
        other_options_samtools_sort = "",
        extra = r"'@RG\tID:{sample}\tSM:{sample}\tCN:SC\tPL:ILLUMINA'", 
    input:
        fwdreads = rules.trimmomatic.output.forward_reads, 
        revreads = rules.trimmomatic.output.reverse_reads
    output:
        bam = "results/02_Mapping/{sample}_bwa_sorted.bam",
    log:
        output = "results/11_Reports/bwa/{sample}.o",
        error = "results/11_Reports/bwa/{sample}.e"
    shell:
        config["MODULES"]["BWA"]+"\n"+config["MODULES"]["SAMTOOLS"]+"""
            (bwa mem -M -T 0 -t {resources.cpus_per_task} -v 1 -R {params.extra} {params.ref} {input.fwdreads} {input.revreads} |
            samtools view -@ {resources.cpus_per_task} {params.other_options_samtools_view} |
            samtools sort -@ {resources.cpus_per_task} {params.other_options_samtools_sort} -o {output.bam} ) 1> {log.output} 2> {log.error}
        """

###############################################################################
rule mark_duplicates:
    message:
        "Picard MarkDuplicates remove PCR duplicates"
    resources: 
        partition='long',
        cpus_per_task=16,
        mem_mb=16000,
        runtime=600,
    input:
        rules.bwa_mapping.output.bam,
    output:
        bam = "results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam",
        metrics="results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt",
    log:
        "results/11_Reports/markduplicatesspark/{sample}_bwa_sorted-mark-dup.log",
    params:
        other_options = "-CREATE_INDEX TRUE -VALIDATION_STRINGENCY SILENT"
    shell:
        config["MODULES"]["PICARDTOOLS"]+"""
            picard MarkDuplicates {params.other_options} -I {input} -O {output.bam} -M {output.metrics} > {log} 2>&1
        """

###############################################################################
rule samtools_index_markdup:
    message:
        "SamTools indexing marked as duplicate BAM file"
    resources: 
        partition='fast',
        cpus_per_task=4,
        mem_mb=4000,
        runtime=120,
    input:
        markdup = rules.mark_duplicates.output.bam
    output:
        index   = "results/02_Mapping/{sample}_bwa_sorted-mark-dup.bai",
    log:
        "results/11_Reports/samtools/{sample}_bwa_sorted-mark-dup-index.log",
    shell:
        config["MODULES"]["SAMTOOLS"]+"""
            samtools index -@ {resources.cpus_per_task} -b {input.markdup} {output.index} &> {log}
        """