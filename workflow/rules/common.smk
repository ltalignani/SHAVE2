#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 common.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile common.smk --cores X --use-conda
# Latest modification:  2023.07.10
# Done:                 
###############################################################################
# CONFIGURATION FILES #
configfile: "config/config.yaml"
cluster_config: "slurm/config.yaml"

import os, sys
from snakemake.utils import min_version
min_version("6.12.0")

from pathlib import Path

###############################################################################
# WILDCARDS #
SAMPLE, = glob_wildcards("raw/{sample}_R1.fastq.gz")
CHROM = "2L 2R 3L 3R X".split()

###############################################################################
# RESOURCES #
TMPDIR = config["resources"]["tmpdir"] # Temporary directory
reference_file =  config["DATA"]["directories"]["reference_file"] #"resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
basename_reference = Path(reference_file).stem

###############################################################################
# FUNCTIONS #

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

###### Config file and sample sheets #####
configfile: "config/config.yaml"

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(samples["sample_name"]), # output: 'A|B|C|D|E'
    unit="|".join(units["unit_name"]), # output: 'lane1|lane2|lane1|lane1|lane1|lane1'


##### Helper functions #####
def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired

def get_fq(wildcards):
    if config["trimming"]["activate"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "results/trimmomatic/{sample}_{group}.fastq.gz",
                        group=["R1", "R2"],
                        **wildcards,
                    ),
                )
            )
        # single end sample
        return {"fq1": "results/trimmomatic/{sample}_single.fastq.gz".format(**wildcards)}
    else:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit)] # u = units.loc[('A', 'lane1')] renvoie toute la ligne de la table unit.tsv:
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{u.fq1}"}
        else:
            return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}



# contigs in reference genome
# def get_contigs():
#     with checkpoints.genome_faidx.get().output[0].open() as fai:
#         return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)


# def get_fastq(wildcards):
#     """Get fastq files of given sample-unit."""
#     fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
#     if len(fastqs) == 2:
#         return {"r1": fastqs.fq1, "r2": fastqs.fq2}
#     return {"r1": fastqs.fq1}


# def is_single_end(sample, unit):
#     """Return True if sample-unit is single end."""
#     return pd.isnull(units.loc[(sample, unit), "fq2"])


# def get_read_group(wildcards):
#     """Denote sample name and platform in read group."""
#     return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
#         sample=wildcards.sample,
#         platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
#     )


# def get_trimmed_reads(wildcards):
#     """Get trimmed reads of given sample-unit."""
#     if not is_single_end(**wildcards):
#         # paired-end sample
#         return expand(
#             "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
#             group=[1, 2],
#             **wildcards
#         )
#     # single end sample
#     return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


# def get_sample_bams(wildcards):
#     """Get all aligned reads of given sample."""
#     return expand(
#         "results/recal/{sample}-{unit}.bam",
#         sample=wildcards.sample,
#         unit=units.loc[wildcards.sample].unit,
#     )


# def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
#     if regions:
#         params = "--intervals '{}' ".format(regions)
#         padding = config["processing"].get("region-padding")
#         if padding:
#             params += "--interval-padding {}".format(padding)
#         return params
#     return default


# def get_call_variants_params(wildcards, input):
#     return (
#         get_regions_param(
#             regions=input.regions, default="--intervals {}".format(wildcards.contig)
#         )
#         + config["params"]["gatk"]["HaplotypeCaller"]
#     )


# def get_recal_input(bai=False):
#     # case 1: no duplicate removal
#     f = "results/mapped/{sample}-{unit}.sorted.bam"
#     if config["processing"]["remove-duplicates"]:
#         # case 2: remove duplicates
#         f = "results/dedup/{sample}-{unit}.bam"
#     if bai:
#         if config["processing"].get("restrict-regions"):
#             # case 3: need an index because random access is required
#             f += ".bai"
#             return f
#         else:
#             # case 4: no index needed
#             return []
#     else:
#         return f


# def get_snpeff_reference():
#     return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


# def get_vartype_arg(wildcards):
#     return "--select-type-to-include {}".format(
#         "SNP" if wildcards.vartype == "snvs" else "INDEL"
#     )


# def get_filter(wildcards):
#     return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}