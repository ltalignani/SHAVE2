######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave2.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile shave2.smk --cores X --use-conda
# Latest modification:  2023.07.10
# Done:                 Updated rules
#
###############################################################################
# PUBLICATIONS #

###############################################################################
# CONFIGURATION #
configfile: "config/config.yaml"
cluster_config = "slurm/config.yaml"

from pathlib import Path
import pandas as pd
from snakemake.utils import min_version
min_version("6.12.0")

import shutil

###############################################################################
# WILDCARDS #
# SAMPLE, READS, = glob_wildcards('raw/{sample}_R{read}.fastq')
SAMPLE, = glob_wildcards("raw/{sample}_R1.fastq.gz"),
CHROM = "2L 2R 3L 3R X".split()
READS = ['1', '2']

###############################################################################
# RESOURCES #
TMPDIR = config["resources"]["tmpdir"] # Temporary directory

###############################################################################
# PARAMETERS #
CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes

reference_file =  config["ref"]["reference"]
basename_reference = Path(reference_file).stem # Pathlib PurePath.stem: remove suffix

###############################################################################
# FUNCTIONS AND COMMANDS #

################################ O N S T A R T #################################
onstart:
    shell("mkdir -p Cluster_logs/")

################################ M O D U L E ####################################
include: "rules/mapping.smk"
include: "rules/polishing.smk"
include: "rules/calling.smk"
include: "rules/vcf_stats.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"

###################### R U L E   D E C L A R A T I O N #########################

################################## A L L #######################################
rule all:
    input:
        fastqc          =   expand(os.path.join("results/00_Quality_Control/fastqc/", '{sample}_R{read}_fastqc.html'), sample=SAMPLE, read=READS),
        fastqscreen     =   "results/00_Quality_Control/fastq-screen/",
        forward_reads   =   expand("results/01_Trimmimg/{sample}_trimmomatic_R1.fastq.gz", sample=SAMPLE),
        reverse_reads   =   expand("results/01_Trimmimg/{sample}_trimmomatic_R2.fastq.gz", sample=SAMPLE),
        forwardUnpaired =   expand("results/01_Trimmimg/{sample}_trimmomatic_unpaired_R1.fastq.gz", sample=SAMPLE),
        reverseUnpaired =   expand("results/01_Trimmimg/{sample}_trimmomatic_unpaired_R2.fastq.gz", sample=SAMPLE),
        bam             = expand("results/02_Mapping/{sample}_bwa_sorted.bam", sample=SAMPLE),
        markdup         = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bam", sample=SAMPLE),
        metrics         = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup_metrics.txt", sample=SAMPLE),
        index_md        = expand("results/02_Mapping/{sample}_bwa_sorted-mark-dup.bai", sample=SAMPLE),
        fixed           = expand("results/04_Polishing/{sample}_md_fixed.bam", sample=SAMPLE),
        index_fx        = expand("results/04_Polishing/{sample}_md_fixed.bai", sample=SAMPLE),
        gvcf            = expand("results/05_Variants/{sample}.{chromosomes}.g.vcf", sample=SAMPLE, chromosomes=CHROM),
        vcf_report      = expand("results/report_vcf.{chromosomes}.html", chromosomes=CHROM),
        db              = expand("results/05_Variants/DB_{chromosomes}", chromosomes=CHROM),
        combined        = expand("results/05_Variants/{chromosomes}_combinedGVCF.vcf.gz", chromosomes=CHROM),
        vcf_file        = expand("results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs.vcf", chromosomes=CHROM),
        vcf_file_combined   = expand("results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs_combined.vcf", chromosomes=CHROM),
        vcf_gz          = expand("results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs.vcf.gz", chromosomes=CHROM),
        tbi_gz          = expand("results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs.vcf.gz.tbi", chromosomes=CHROM),
        vcf_filt        = expand("results/05_Variants/filtered/All_samples.{chromosomes}.GenotypeGVCFs.hf.vcf.gz", chromosomes=CHROM),
        check           = expand("results/00_Quality_Control/validatesamfile/{sample}_md_fixed_ValidateSam.txt", sample=SAMPLE),
        flagstat        = expand("results/00_Quality_Control/{sample}_md_fixed_bam.flagstat.txt", sample=SAMPLE),
        idxstats        = expand("results/00_Quality_Control/{sample}_md_fixed.idxstats.txt", sample=SAMPLE),
        stats           = expand("results/00_Quality_Control/stats/{sample}_md_fixed_stats.txt", sample=SAMPLE),


