#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 filtering.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile shave2.smk --cores X --use-conda
# Latest modification:  2023.02.02
# Done:                 Divided snakefile rules
###############################################################################
rule gatk_filter:
    # Aim: Filter variant calls based on INFO and/or FORMAT annotations.
    message: "Hard-filtering of vcf files"
    resources: 
        partition='long',
        mem_mb=12000,
        runtime=600,
    input:
        ref=reference_file,
        vcf_file_all =rules.bcftools_concat.output.vcf_gz,
    output:
        vcf_file = "results/05_Variants/filtered/All_samples.{chromosomes}.GenotypeGVCFs.hf.vcf.gz",
    params: 
        filters={"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="",

    log:
        "results/11_Reports/variantfiltration/{chromosomes}.merged_hardfiltered.log",
    shell:
        config["MODULES"]["GATK4"]+"""
            gatk VariantFiltration -R {input.ref} -V {input.vcf_file_all} -O {output.vcf_file} --filterExpression {params.filters} --filterName 'my_filters'
        """
