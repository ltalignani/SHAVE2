#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 genotypegvcfs.smk
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
from pathlib import Path
import pandas as pd
from snakemake.utils import min_version
min_version("6.12.0")

import shutil

###############################################################################
# WILDCARDS #
# SAMPLE, = glob_wildcards("raw/{sample}_R1.fastq.gz"),
# CHROM = "2L 2R 3L 3R X".split()
# READS = ['1', '2']

###############################################################################
# RESOURCES #
TMPDIR = config["resources"]["tmpdir"] # Temporary directory,

###############################################################################
# PARAMETERS #
CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes

###############################################################################
# FUNCTIONS #

###############################################################################
rule create_sequence_dict:
    message: "create sequence dict for gatk_HaplotypeCaller reference"
    resources: 
        partition='fast',
        mem_mb=8000,
        runtime=120,
    input:
        reference = reference_file
    output:
        dictionary = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.dict"
    log:
        error =  f'results/11_Reports/create_sequence_dict/{basename_reference}.e',
        output = f'results/11_Reports/create_sequence_dict/{basename_reference}.o'
    shell: config["MODULES"]["GATK4"]+"""
        gatk CreateSequenceDictionary --java-options "-Xmx{resources.mem_mb}M" -R {input.reference} -O {output.dictionary} 1>{log.output} 2>{log.error}
    """

###############################################################################
rule create_sequence_faidx:
    message: "create sequence fai for gatk_HaplotypeCaller reference"
    resources: 
        partition='fast',
        mem_mb=8000,
        runtime=120,
    input:
        reference = reference_file
    output:
        fai = "resources/genomes/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fasta.fai"
    log:
        error =  f'results/11_Reports/samtools/create_sequence_fai/{basename_reference}.e',
        output = f'results/11_Reports/samtools/create_sequence_fai/{basename_reference}.o'

    shell: config["MODULES"]["SAMTOOLS"]+"""
        samtools faidx {input.reference} 1>{log.output} 2>{log.error}
    """

###############################################################################
rule HaplotypeCaller:
    message:
        "GATK's HaplotypeCaller SNPs and indels calling for {wildcards.sample} sample"
    resources: 
        partition='long',
        mem_mb=20000,
        runtime=24000,
    input:
        bam = rules.fixmateinformation.output.fixed,
        reference = reference_file,
        dictionary = rules.create_sequence_dict.output.dictionary,
        fai = rules.create_sequence_faidx.output.fai,
    output:
        gvcf = "results/05_Variants/{sample}.{chromosomes}.g.vcf",
   
    log:
        output = "results/11_Reports/haplotypecaller/{sample}.{chromosomes}_variant-call.o",
        error = "results/11_Reports/haplotypecaller/{sample}.{chromosomes}_variant-call.e",
    params: 
        other_options = config["gatk"]["haplotypecaller"], # -ERC GVCF
        # alleles = config["alleles"][alleles_target] 
        interval = "{chromosomes}",
    benchmark:
        "benchmarks/haplotypecaller/{sample}_{chromosomes}_variant-call.tsv"
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk HaplotypeCaller --java-options "{params.java_opts} -Xmx{resources.mem_mb}M" -nct 1 -R {input.reference} -I {input.bam} -O {output.gvcf} \
            {params.other_options} \
            --native-pair-hmm-threads 1 \
            -L {params.interval} \
            1>{log.output} 2>{log.error}
        """

###############################################################################
def get_gvcf_list(list):
    gvcf_list = " -V ".join(list)
    return(f"-V {gvcf_list}")

##############################################################################
rule GenomicsDBImport:
    message:
        "GATK's GenomicsDBImport for multiple g.vcfs for chromosome {wildcards.chromosomes}"
    resources: 
        partition='long',
        mem_mb=10000,
        runtime=24000,  
    input:
        gvcf_list   = expand(rules.HaplotypeCaller.output.gvcf, sample=SAMPLE, chromosomes=CHROM),
        reference   = reference_file,
    output:
        db          = directory("results/05_Variants/DB_{chromosomes}"),
        combined    = "results/05_Variants/{chromosomes}_combinedGVCF.vcf.gz"
    log:
        output      = "results/11_Reports/genomicsdbimport/{chromosomes}_genomicsdbimport.o",
        error       = "results/11_Reports/genomicsdbimport/{chromosomes}_genomicsdbimport.e", 
    params:
        interval    = "{chromosomes}",
        str_join    = get_gvcf_list(expand(rules.HaplotypeCaller.output.gvcf, sample=SAMPLE, chromosomes=CHROM)),
        other_options = config["gatk"]["genomicsdbimport"]
    shell:
        config["MODULES"]["GATK4"]+"""
        gatk GenomicsDBImport --java-options "-Xmx{resources.mem_mb}M" -R {input.reference} {params.str_join} \
            --genomicsdb-workspace-path {output.db} \ 
            {params.other_options} \
            -L {params.interval} \
            1>{log.output} 2>{log.error}

        gatk CombineGVCFs --java-options "-Xmx{resources.mem_mb}M" -R {input.reference} {params.str_join} \
        -L {params.interval} -O {output.combined} 1>>{log.output} 2>>{log.error}
        """

##############################################################################
rule GenotypeGVCFs_merge:
    message: "GATK's GenotypeGVCFs for chromosome {wildcards.chromosomes}"
    resources:
        partition='long',
        mem_mb=40000, 
        runtime=24000,    
    input:
        genomicsdb          = rules.GenomicsDBImport.output.db,
        combined            = rules.GenomicsDBImport.output.combined,
        reference           = reference_file,
    output:
        vcf_file            = "results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs.vcf",
        vcf_file_combined   = "results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs_combined.vcf",
    params:
        interval            = "{chromosomes}",
        other_options       = config["gatk"]["genotypegvcfs"],
    log:
        error               = "results/11_Reports/genotypegvcfs/gatk_GenotypeGVCFs.{chromosomes}.merge/bwa.e",
        output              = "results/11_Reports/genotypegvcfs/gatk_GenotypeGVCFs.{chromosomes}.merge/bwa.o"
    benchmark:
        "benchmarks/genotypegvcfs/{chromosomes}_bwa_genotyped.tsv"
    shell:
        config["MODULES"]["GATK4"]+"""
        gatk GenotypeGVCFs --java-options "-Xmx{resources.mem_mb}M" -R {input.reference} -V gendb://{input.genomicsdb} -L {params.interval} -O {output.vcf_file} {params.other_options} 1>{log.output} 2>{log.error}
        gatk GenotypeGVCFs --java-options "-Xmx{resources.mem_mb}M" -R {input.reference} -V {input.combined} -L {params.interval} -O {output.vcf_file_combined} {params.other_options} 1>>{log.output} 2>>{log.error}
        """

##############################################################################
rule bcftools_concat:
    message: "Concatenate vcfs produced for each interval"
    resources: 
        partition='fast',
        cpus_per_task=8,
        mem_mb=20000,
        runtime=12000, 
    input:
        vcf_file_all = rules.GenotypeGVCFs_merge.output.vcf_file,
    output:
        vcf_gz = "results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs.vcf.gz",
        tbi_gz = "results/05_Variants/genotypegvcfs/All_samples.{chromosomes}.GenotypeGVCFs.vcf.gz.tbi",
    log:
        error =  "results/11_Reports/bcftools_concat/bcftools_concat.{chromosomes}.e",
        output = "results/11_Reports/bcftools_concat/bcftools_concat.{chromosomes}.o"
    shell:
        config["MODULES"]["BCFTOOLS"]+"""
        bcftools concat --threads {resources.cpus_per_task} {input.vcf_file_all} -Oz -o {output.vcf_gz} 1>{log.output} 2>{log.error}
        bcftools index --threads {resources.cpus_per_task} --tbi {output.vcf_file} 1>>{log.output} 2>>{log.error}
        """
