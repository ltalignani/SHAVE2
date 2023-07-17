#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 vcf_stats.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile vcf_stats.smk --cores X --use-conda
# Latest modification:  2023.02.02
# Done:                 Divided snakefile rules
###############################################################################
rule vcf_stats:
    message: "Execute rule VCF stats"
    input:
        vcf_file_all = rules.bcftools_concat.output.vcf_gz
    output:
        freq = 'results/06_SNP_calling_stats/All_samples.{chromosomes}.GenotypeGVCFs.frq',
        depth = 'results/06_SNP_calling_stats/All_samples.{chromosomes}.GenotypeGVCFs.idepth',
        depth_mean = 'results/06_SNP_calling_stats/All_samples.{chromosomes}.GenotypeGVCFs.ldepth.mean',
        qual = 'results/06_SNP_calling_stats/All_samples.{chromosomes}.GenotypeGVCFs.lqual',
        missing_ind = 'results/06_SNP_calling_stats/All_samples.{chromosomes}.GenotypeGVCFs.imiss',
        miss = 'results/06_SNP_calling_stats/All_samples.{chromosomes}.GenotypeGVCFs.lmiss',
    resources: cpus=1, mem_mb=4000, time_min=120
    params: partition = 'long',
    log:
            error =  'results/11_Reports/vcf_stats/vcftools.{chromosomes}.e',
            output = 'results/11_Reports/vcf_stats/vcftools.{chromosomes}.o'        
    shell:
        config["MODULES"]["VCFTOOLS"]+" && "+config["MODULES"]["SAMTOOLS"]+"""
            vcftools --gzvcf {input.vcf_file_all}  --remove-indels --freq2 --max-alleles 3 --stdout 1> {output.freq}
            vcftools --gzvcf {input.vcf_file_all}  --remove-indels --depth --stdout 1> {output.depth}
            vcftools --gzvcf {input.vcf_file_all}  --remove-indels --site-mean-depth --stdout 1> {output.depth_mean}
            vcftools --gzvcf {input.vcf_file_all}  --remove-indels --site-quality --stdout 1> {output.qual}
            vcftools --gzvcf {input.vcf_file_all}  --remove-indels --missing-indv --stdout 1> {output.missing_ind}
            vcftools --gzvcf {input.vcf_file_all}  --remove-indels --missing-site --stdout 1> {output.miss}
        """
#cat {output.missing_ind} | awk \'\{if($5>0.75)print $1\}\'|grep -v INDV> remove-indv_75perc

###############################################################################
rule report_vcf:
    message: "Execute rule report_vcf",
    resources: cpus=1, mem_mb=4000, time_min=1200
    input:
        freq = rules.vcf_stats.output.freq,
        depth = rules.vcf_stats.output.depth,
        depth_mean = rules.vcf_stats.output.depth_mean,
        qual = rules.vcf_stats.output.qual,
        missing_ind = rules.vcf_stats.output.missing_ind,
        miss = rules.vcf_stats.output.miss
    params:
        outdir = 'results/',
        partition = 'long',
    output:
        report = "results/report_vcf.{chromosomes}.html",
    log:
        error =  'results/11_Reports/report_vcf/report_vcf.{chromosomes}.e',
        output = 'results/11_Reports/report_vcf/report_vcf.{chromosomes}.o'
    script:
        """workflow/scripts/report_vcf.Rmd"""

###############################################################################