#######M######I#####V#####E######G######E######C######I#######R######D#########
# Name:                 shave2.smk
# Author:               Loïc TALIGNANI
# Affiliation:          IRD_UMR224_MIVEGEC
# Aim:                  Snakefile for SHort-read Alignment for VEctor pipeline
# Date:                 2022.10.05
# Run:                  snakemake --snakefile shave2.smk --cores X --use-conda
# Latest modification:  2022.11.28
# Done:                 Added HaplotypeCaller, GenotypeGVCFs, VariantFiltration
#                       java Garbage Collection limit
###############################################################################
# PUBLICATIONS #

###############################################################################
# CONFIGURATION #
configfile: "config/config.yaml"

from snakemake.utils import min_version
min_version("5.18.0")

import shutil

###############################################################################
# WILDCARDS #
SAMPLE, = glob_wildcards("resources/reads/{sample}_R1.fastq.gz")
sample_ids, aligner_ids, markdup_ids, mincov_ids = glob_wildcards("results/04_Variants/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.g.vcf")

###############################################################################
# RESOURCES #
OS = config["os"]                      # Operating system
CPUS = config["resources"]["cpus"]     # Threads
MEM_GB = config["resources"]["mem_gb"] # Memory (RAM) in Gb
TMPDIR = config["resources"]["tmpdir"] # Temporary directory

###############################################################################
# ENVIRONMENTS #
FASTQC = config["conda"][OS]["fastqc"]              # FastQC
FASTQSCREEN = config["conda"][OS]["fastq-screen"]   # Fastq-Screen
MULTIQC = config["conda"][OS]["multiqc"]            # MultiQC
CUTADAPT = config["conda"][OS]["cutadapt"]          # Cutadapt
SICKLETRIM = config["conda"][OS]["sickle-trim"]     # Sickle-trim
BOWTIE2 = config["conda"][OS]["bowtie2"]            # Bowtie2
BWA = config["conda"][OS]["bwa"]                    # Bwa
SAMTOOLS = config["conda"][OS]["samtools"]          # SamTools
BEDTOOLS = config["conda"][OS]["bedtools"]          # BedTools
BCFTOOLS = config["conda"][OS]["bcftools"]          # BcfTools
GAWK = config["conda"][OS]["gawk"]                  # Gawk
LOFREQ = config["conda"][OS]["lofreq"]              # LoFreq
GATK = config["conda"][OS]["gatk"]                  # GATK 3.8
GATK4 = config["conda"][OS]["gatk4"]                # GATK 4.3.0
PICARD = config["conda"][OS]["picard"]              # Picard 2.24.7

###############################################################################
# PARAMETERS #
LENGTHc = config["cutadapt"]["length"]              # Cutadapt --minimum-length
TRUSEQ = config["cutadapt"]["kits"]["truseq"]       # Cutadapt --adapter Illumina TruSeq
NEXTERA = config["cutadapt"]["kits"]["nextera"]     # Cutadapt --adapter Illumina Nextera
SMALL = config["cutadapt"]["kits"]["small"]         # Cutadapt --adapter Illumina Small

COMMAND = config["sickle-trim"]["command"]          # Sickle-trim command
ENCODING = config["sickle-trim"]["encoding"]        # Sickle-trim --qual-type
QUALITY = config["sickle-trim"]["quality"]          # Sickle-trim --qual-threshold
LENGTH = config["sickle-trim"]["length"]            # Sickle-trim --length-treshold

CONFIG = config["fastq-screen"]["config"]           # Fastq-screen --conf
MAPPER = config["fastq-screen"]["aligner"]          # Fastq-screen --aligner
SUBSET = config["fastq-screen"]["subset"]           # Fastq-screen --subset

ALIGNER = config["aligner"]                         # Aligners ('bwa' or 'bowtie2')
MARKDUP = config["markdup"]                         # Mark Duplicate Program ('picard' or 'samtools')

BWAPATH = config["bwa"]["path"]                     # BWA path to indexes
BT2PATH = config["bowtie2"]["path"]                 # Bowtie2 path to indexes
SENSITIVITY = config["bowtie2"]["sensitivity"]      # Bowtie2 sensitivity preset

REFPATH = config["consensus"]["path"]               # Path to genome reference
REFERENCE = config["consensus"]["reference"]        # Genome reference sequence, in fasta format
MINCOV = config["consensus"]["mincov"]              # Minimum coverage, mask lower regions with 'N'
MINAF = config["consensus"]["minaf"]                # Minimum allele frequency allowed
IUPAC = config["consensus"]["iupac"]                # Output variants in the form of IUPAC ambiguity codes

###############################################################################
# FUNCTIONS AND COMMANDS#

############################### O N S T A R T ################################
onstart:
    print("##### Creating profile pipeline #####\n") 
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p Cluster_logs/haplotypecaller")
    shell("mkdir -p Cluster_logs/bwa_mapping")

############################# O N S U C C E S S ##############################
onsuccess:
    shutil.rmtree(".snakemake")

################################## A L L #####################################
rule all:
    input:
        multiqc = "results/00_Quality_Control/multiqc/",
        fastqc = "results/00_Quality_Control/fastqc/",
        fastqscreen = "results/00_Quality_Control/fastq-screen/",
        covstats = expand("results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_coverage-stats.tsv",
                          sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        index_archive = expand("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_variant-filt.gz.tbi",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),        
        archive = expand("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_variant-filt.vcf.gz",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        vcfarchive = expand("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.vcf.gz",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),       
        index = expand("results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bai",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),        
        stats = expand("results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.txt",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        callable_loci = expand("results/05_Validation/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate_callable_status.bed",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV),
        check = expand("results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam",
                           sample = SAMPLE, aligner = ALIGNER, markdup = MARKDUP, mincov = MINCOV)

###############################################################################
rule tabix_tabarch_indexing:
    # Aim: tab archive indexing
    # Use: tabix [OPTIONS] [TAB.bgz]
    message:
        "Tabix tab archive indexing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        SAMTOOLS
    input:
        archive = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_variant-filt.vcf.gz"
    output:
        index = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_variant-filt.gz.tbi"
    log:
        "results/11_Reports/tabix/{sample}_{aligner}_{markdup}_{mincov}X_variant-archive-index.log"
    shell:
        "tabix "             # Tabix, indexes a TAB-delimited genome position file in.tab.bgz and creates an index file
        "{input.archive} "   # The input data file must be position sorted and compressed by bgzip
        "-f "                # overwrite index if already existing   
        "1> {output.index} " # Tabix output TBI index formats
        "2> {log}"           # Log redirection

###############################################################################
rule bcftools_variant_filt_archive:
    # Aim: Variant block compressing
    # Use: bgzip [OPTIONS] -c -@ [THREADS] [INDEL.vcf] 1> [COMPRESS.vcf.bgz]
    message:
        "Bgzip variant block compressing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BCFTOOLS
    threads:
        CPUS
    input:
        variantfilt = "results/04_Variants/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.vcf"                  #results/04_Variants/lofreq/{sample}_{aligner}_{mincov}X_variant-filt.vcf"
    output:
        archive = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_variant-filt.vcf.gz"
    log:
        "results/11_Reports/bcftools/{sample}_{aligner}_{markdup}_{mincov}X_variant-archive.log"
    shell:
        "bcftools "                         # bcftools,  a set of utilities that manipulate variant calls in the Variant Call Format (VCF).
        "view "                             # view : subset, filter and convert VCF and BCF files
        "--threads {threads} "              # -@: Number of threads to use (default: 1)
        "{input.variantfilt} "              # VCF input file,
        "-Oz -o {output.archive} "          # -O[z|b]: output-type -o: VCF output file,
        "&> {log}"                          # Log redirection

###############################################################################
rule hard_filter_calls:
    # Aim: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller.
    # In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
    # Use: gatk --java-options "-Xmx4g" GenotypeGVCFs \
    #      -R Homo_sapiens_assembly38.fasta \
    #      -V input.g.vcf.gz \
    #      -O output.vcf.gz
    message:
        "Hard-filtering for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GATK4
    input:
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        vcf="results/04_Variants/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.vcf",
    output:
        vcf=temp("results/04_Variants/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.vcf"),
    params:
        filters={"myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"},
        extra="",
        java_opts="-XX:ParallelGCThreads=8 -XX:ConcGCThreads=2 -XX:G1ConcRefinementThreads=9",
    resources:
        mem_gb= MEM_GB
    threads: 
        CPUS
    benchmark:
        "benchmarks/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.tsv"
    log:
        "results/11_Reports/variantfiltration/{sample}_{aligner}_{markdup}_{mincov}X_hardfiltered.log",
    wrapper:
        "0.74.0/bio/gatk/variantfiltration"

###############################################################################
rule bcftools_genotype_gvcfs_archive:
    # Aim: Variant block compressing
    # Use: bcftools view -@ [THREADS] [input.vcf] -Oz -o [COMPRESS.vcf.gz]
    message:
        "Bgzip variant block compressing for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BCFTOOLS
    threads:
        CPUS
    input:
        vcf = "results/04_Variants/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.vcf"                  #results/04_Variants/lofreq/{sample}_{aligner}_{mincov}X_variant-filt.vcf"
    output:
        vcfarchive = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.vcf.gz"
    log:
        "results/11_Reports/bcftools/{sample}_{aligner}_{markdup}_{mincov}X_variant-archive.log"
    shell:
        "bcftools "                         # bcftools,  a set of utilities that manipulate variant calls in the Variant Call Format (VCF).
        "view "                             # view : subset, filter and convert VCF and BCF files
        "--threads {threads} "              # -@: Number of threads to use (default: 1)
        "{input.vcf} "                      # VCF input file, 
        "-Oz -o {output.vcfarchive} "       # -O[z|b]: output-type -o: VCF output file,
        "&> {log}"                          # Log redirection

###############################################################################
rule genotype_single_gvcfs:
    # Aim: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller.
    # In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
    # Use: gatk --java-options "-Xmx4g" GenotypeGVCFs \
    #      -R Homo_sapiens_assembly38.fasta \
    #      -V input.g.vcf.gz \
    #      -O output.vcf.gz
    input:
        gvcf="results/04_Variants/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.g.vcf",
        # N.B. gvcf or genomicsdb must be specified
        # in the latter case, this is a GenomicsDB data store
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    output:
        vcf=temp("results/04_Variants/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.vcf")
    log:
        "results/11_Reports/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotypegvcfs.log"
    params:
        extra="",  # optional
        java_opts="-XX:ParallelGCThreads=8 -XX:ConcGCThreads=2 -XX:G1ConcRefinementThreads=9", # optional
    resources:
        mem_mb=16000
    benchmark:
        "benchmarks/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.tsv"
    wrapper:
        "v1.16.0/bio/gatk/genotypegvcfs"

###############################################################################
rule genotype_multiple_gvcfs:
    # Aim: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller.
    # In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
    # Use: gatk --java-options "-Xmx4g" GenotypeGVCFs \
    #      -R Homo_sapiens_assembly38.fasta \
    #      -V input.g.vcf.gz \
    #      -O output.vcf.gz
    input:
        gvcf= directory("results/04_Variants/genomics_db"), # combined gvcf over multiple samples "results/04_Variants/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.g.vcf",
        # N.B. gvcf or genomicsdb must be specified
        # in the latter case, this is a GenomicsDB data store
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    output:
        vcf=temp("results/04_Variants/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.vcf")
    log:
        "results/11_Reports/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotypegvcfs.log"
    params:
        extra="",  # optional
        java_opts="-XX:ParallelGCThreads=8 -XX:ConcGCThreads=2 -XX:G1ConcRefinementThreads=9", # optional
    resources:
        mem_mb=16000
    benchmark:
        "benchmarks/genotypegvcfs/{sample}_{aligner}_{markdup}_{mincov}X_genotyped.tsv"
    wrapper:
        "v1.16.0/bio/gatk/genotypegvcfs"

###############################################################################
rule genomics_db_import:
    # Aim: Import single-sample GVCFs into GenomicsDB before joint genotyping.
    # The GATK4 Best Practice Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from multiple samples. GenomicsDBImport 
    #  offers the same functionality as CombineGVCFs and comes from the Intel-Broad Center for Genomics. The datastore transposes 
    # sample-centric variant information across genomic loci to make data more accessible to tools.
    # Use: gatk --java-options "-Xmx4g -Xms4g" \
    #      GenomicsDBImport \
    #      --genomicsdb-workspace-path my_database \
    #      --batch-size 50 \
    #      -L chr1:1000-10000 \
    #      --sample-name-map cohort.sample_map \
    #      --tmp-dir=/path/to/large/tmp \
    #      --reader-threads 5
    message:
        "GATK's GenomicsDBImport for multiple g.vcfs"
    input:
        expand("results/04_Variants/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.g.vcf",
                sample=sample_ids, aligner=aligner_ids, markdup=markdup_ids, mincov=mincov_ids)        #["calls/a.g.vcf.gz", "calls/b.g.vcf.gz"], #"results/04_Variants/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.g.vcf"
    output:
        db=directory("results/04_Variants/genomics_db"),
    log:
        "results/11_Reports/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_genomicsdbimport.log",
    params:
        intervals="ref",
        db_action="create",  # optional
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=16000,
    wrapper:
        "v1.19.2/bio/gatk/genomicsdbimport"

###############################################################################
rule haplotype_caller_gvcf:
    # Aim: Call germline SNPs and indels via local re-assembly of haplotypes
    # Use: gatk --java-options "-Xmx4g" HaplotypeCaller  \
    #      -R Homo_sapiens_assembly38.fasta \
    #      -I input.bam \
    #      -O output.g.vcf.gz \
    #      -ERC GVCF or --emit-ref-confidence GVCF
    message:
        "GATK's HaplotypeCaller SNPs and indels calling for {wildcards.sample} sample ({wildcards.aligner})"
    input:
        # single or list of bam files
        bam="results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bam",
        index = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bai",
        ref="resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        # known="dbsnp.vcf"  # optional
    output:
        gvcf=temp("results/04_Variants/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.g.vcf"),
    #       bam="{sample}.assemb_haplo.bam",
    log:
        "results/11_Reports/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.log",
    params:
        extra="", 
        java_opts="-XX:ParallelGCThreads=8 -XX:ConcGCThreads=2 -XX:G1ConcRefinementThreads=9" ,  # -XX:ParallelGCThreads=NB_CPUS -XX:ConcGCThreads=NB_CPUS/4 -XX:G1ConcRefinementThreads=NB_CPUS+1
    threads: CPUS
    resources:
        mem_mb=16000,
    benchmark:
        "benchmarks/haplotypecaller/{sample}_{aligner}_{markdup}_{mincov}X_variant-call.tsv"
    wrapper:
        "v1.16.0/bio/gatk/haplotypecaller"

###############################################################################
rule samtools_stats:
    # Aim: Collects statistics from BAM files
    # Use: samtools stats -r ref.fa input.bam
    message:
        "SamTools indexing marked as duplicate BAM file {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
        CPUS
    input:
        bam = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bam",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta"
    output:
        stats = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.txt"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.log"
    shell:
        "samtools stats "                                                   # Samtools stats, collects statistics from BAM files. The output can be visualized using plot-bamstats.
        "--threads {threads} "                                              # -@: Number of additional threads to use (default: 1)
        "-r {input.refpath} "                                               # -r: Reference sequence (required for GC-depth and mismatches-per-cycle calculation).
        "{input.bam} "                                                      # mark-dup bam input
        "1> {output.stats} "                                                # stats output
        "2> {log}"                                                          # Log redirection

###############################################################################
rule callable_loci:
    # Aim: Collects statistics on callable, uncallable, poorly mapped, and other parts of the genome.
    # A very common question about a NGS set of reads is what areas of the genome are considered callable. This tool
    # considers the coverage at each locus and emits either a per base state or a summary interval BED file that
    # partitions the genomic intervals into the following callable states:
    # REF_N: The reference base was an N, which is not considered callable the GATK
    # PASS: The base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE
    # NO_COVERAGE: Absolutely no reads were seen at this locus, regardless of the filtering parameters
    # LOW_COVERAGE: There were fewer than min. depth bases at the locus, after applying filters
    # EXCESSIVE_COVERAGE: More than -maxDepth read at the locus, indicating some sort of mapping problem
    # POOR_MAPPING_QUALITY: More than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads
    # Use: gatk3 -T CallableLoci \
    #     -T CallableLoci \
    #     -R reference.fasta \
    #     -I myreads.bam \
    #     -summary table.txt \
    #     -o callable_status.bed
    message:
        "GATK3 CallableLoci for {wildcards.sample} sample ({wildcards.aligner}"
    conda:
        GATK
    input:
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
        bam = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bam",
        index = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bai",
    output:
        call = "results/05_Validation/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate_callable_status.bed",
        summary = "results/05_Validation/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate_summary_table.txt"
    log :
        "results/11_reports/callableloci/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate_callable_status.log"
    shell:
        "gatk3 -T CallableLoci -R {input.refpath} -I {input.bam} -summary {output.summary} -o {output.call}" #  > {log} 2>&1

###############################################################################
rule validate_sam:
    # Aim: Basic check for bam file validity, as interpreted by the Broad Institute.
    # Use: picard.jar ValidateSamFile \
    #      -I input.bam \
    #      - MODE SUMMARY
    message:
        "Picard ValidateSamFile for {wildcards.sample} sample ({wildcards.aligner})"
    input:
        bam = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bam",
        index = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bai",
        refpath = "resources/genomes/GCA_018104305.1_AalbF3_genomic.fasta",
    output:
        check = "results/05_Validation/validatesamfile/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.txt"
    log:
        "results/11_Reports/validatesamfiles/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.log"
    shell:
        """
        picard ValidateSamFile -I {input.bam} -R {input.refpath} -O {output.check} --VERBOSITY ERROR > {log} 2>&1
        """

###############################################################################
rule samtools_indel_indexing:
    # Aim: indexing indel qualities BAM file
    # Use: samtools index -@ [THREADS] -b [INDELQUAL.bam] [INDEX.bai]
    message:
        "SamTools indexing indel qualities BAM file {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        SAMTOOLS
    threads:
       CPUS
    input:
        fixmate = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bam"
    output:
        index = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate-index.log"
    threads: 
        CPUS
    shell:
        "samtools index "      # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {threads} "        # Number of additional threads to use (default: 0)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.fixmate} "   # Sorted bam input
        "{output.index} "      # Markdup bam output
        "&> {log}"             # Log redirection

###############################################################################
rule fixmateinformation:
    # Aim: This tool ensures that all mate-pair information is in sync between each read and its mate pair.
    #      If no #OUTPUT file is supplied then the output is written to a temporary file and then copied over 
    #      the #INPUT file (with the original placed in a .old file.)
    # Use: picard.jar FixMateInformation \
    #      -I input.bam \
    #      -O fixed_mate.bam \
    #      --ADD_MATE_CIGAR true
    message:
        "Picard FixMateInformation for {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        PICARD
    input:
        bam = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam",
    output:
        fixmate = "results/05_Validation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate.bam"
    threads: CPUS
    log:
        "results/11_Reports/fixmateinformation/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual_fixed-mate_bam.log"
    shell:
        """
        picard FixMateInformation -I {input.bam} -O {output.fixmate} --ADD_MATE_CIGAR true
        """

###############################################################################
rule lofreq_indel_qualities:
    # Aim: Indels qualities. Can be used instead of GATK’s BQSR or on non-Illumina data. 
    #      If you have Illumina data and don’t want to use GATK’s BQSR then the easiest thing is to use the --dindel option.
    # Use: lofreq indelqual --dindel -f [MASKEDREF.fasta] -o [INDEL.bam] [MARKDUP.bam]
    # Note: do not realign your BAM file afterwards!
    message:
        "LoFreq insert indels qualities for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        LOFREQ
    input:
        maskedref = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_masked-ref.fasta",
        markdup = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam"
    output:
        indelqual = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.bam"
    log:
        "results/11_Reports/lofreq/{sample}_{aligner}_{markdup}_{mincov}X_indel-qual.log"
    shell:
        "lofreq "                  # LoFreq, fast and sensitive inference of SNVs and Indels
        "indelqual "                # Insert indel qualities into BAM file (required for indel predictions)
        "--dindel "                 # Add Dindel's indel qualities Illumina specifics (need --ref and clashes with -u)
        "--ref {input.maskedref} "  # -f: Reference (masked) sequence used for mapping (only required for --dindel)
        "--out {output.indelqual} " # -o: Indel BAM file output (default: standard output)
        "{input.markdup} "          # Markdup BAM input
        "&> {log}"                  # Log redirection

###############################################################################
rule bedtools_masking:
    # Aim: masking low coverage regions
    # Use: bedtools maskfasta [OPTIONS] -fi [REFERENCE.fasta] -bed [RANGE.bed] -fo [MASKEDREF.fasta]
    message:
        "BedTools masking low coverage regions for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BEDTOOLS
    params:
        path = REFPATH,
        reference = REFERENCE
    input:
        lowcovmask = "results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_low-cov-mask.bed"
    output:
        maskedref = "results/04_Variants/{sample}_{aligner}_{markdup}_{mincov}X_masked-ref.fasta"
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{markdup}_{mincov}X_masking.log"
    shell:
        "bedtools maskfasta "                        # Bedtools maskfasta, mask a fasta file based on feature coordinates
        "-fi {params.path}{params.reference} "       # Input FASTA file
        "-bed {input.lowcovmask} "                   # BED/GFF/VCF file of ranges to mask in -fi
        "-fo {output.maskedref} "                    # Output masked FASTA file
        "&> {log}"                                   # Log redirection

###############################################################################
rule bedtools_merged_mask:
    # Aim: merging overlaps
    # Use: bedtools merge [OPTIONS] -i [FILTERED.bed] -g [GENOME.fasta]
    message:
        "BedTools merging overlaps for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        BEDTOOLS
    input:
        mincovfilt = "results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.bed"
    output:
        lowcovmask = temp("results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_low-cov-mask.bed")
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{markdup}_{mincov}X_merging.log"
    shell:
        "bedtools merge "        # Bedtools merge, merges overlapping BED/GFF/VCF entries into a single interval
        "-i {input.mincovfilt} "  # -i: BED/GFF/VCF input to merge
        "1> {output.lowcovmask} " # merged output
        "2> {log}"                # Log redirection

###############################################################################
rule awk_mincovfilt:
    # Aim: minimum coverage filtration
    # Use: awk '$4 < [MINCOV]' [BEDGRAPH.bed] 1> [FILTERED.bed]
    message:
        "Awk minimum coverage filtration for {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GAWK
    params:
        mincov = MINCOV
    input:
        genomecov = "results/03_Coverage/{sample}_{aligner}_{markdup}-genome-cov.bed"
    output:
        mincovfilt = temp("results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.bed")
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{markdup}_{mincov}X_min-cov-filt.log"
    shell:
        "awk "                      # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "'$4 < {params.mincov}' "    # Minimum coverage for masking regions in consensus sequence
        "{input.genomecov} "         # BedGraph coverage input
        "1> {output.mincovfilt} "    # Minimum coverage filtered bed output
        "2> {log} "                  # Log redirection

###############################################################################
rule awk_coverage_stats:
    # Aim: computing genomme coverage stats
    # Use: awk {FORMULA} END {{print [RESULTS.tsv] [BEDGRAPH.bed]
    message:
        "Awk compute genome coverage statistics BED {wildcards.sample} sample ({wildcards.aligner}-{wildcards.mincov})"
    conda:
        GAWK
    params:
        mincov = MINCOV
    input:
        genomecov = "results/03_Coverage/{sample}_{aligner}_{markdup}-genome-cov.bed"
    output:
        covstats = "results/03_Coverage/{sample}_{aligner}_{markdup}_{mincov}X_coverage-stats.tsv"
    log:
        "results/11_Reports/awk/{sample}_{aligner}_{markdup}_{mincov}X_coverage-stats.log"
    shell:
        "awk ' "                                  # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "$4 >= {params.mincov} "                   # Minimum coverage
        "{{supMinCov+=$3-$2}} ; "                  # Genome size >= @ mincov X
        "{{genomeSize+=$3-$2}} ; "                 # Genome size
        "{{totalBases+=($3-$2)*$4}} ; "            # Total bases @ 1 X
        "{{totalBasesSq+=(($3-$2)*$4)**2}} "       # Total bases square @ 1 X
        "END "                                     # End
        "{{print "                                 # Print
        """ "sample_id", "\t", """                 # Sample ID header
        """ "mean_depth", "\t", """                # Mean depth header
        """ "standard_deviation", "\t", """        # Standard deviation header
        """ "cov_percent_@{wildcards.mincov}X" """ # Coverage percent @ mincov X header
        "ORS "                                     # \n newline
        """ "{wildcards.sample}", "\t", """        # Sample ID value
        """ int(totalBases/genomeSize), "\t", """  # Mean depth value
        """ int(sqrt((totalBasesSq/genomeSize)-(totalBases/genomeSize)**2)), "\t", """ # Standard deviation value
        """ supMinCov/genomeSize*100 """           # Coverage percent @ mincov X value
        "}}' "                                     #
        "{input.genomecov} "                       # BedGraph coverage input
        "1> {output.covstats} "                    # Mean depth output
        "2> {log}"                                 # Log redirection

###############################################################################
rule bedtools_genome_coverage:
    # Aim: computing genome coverage sequencing
    # Use: bedtools genomecov [OPTIONS] -ibam [MARKDUP.bam] 1> [BEDGRAPH.bed]
    message:
        "BedTools computing genome coverage for {wildcards.sample} sample against reference genome sequence ({wildcards.aligner})"
    conda:
        BEDTOOLS
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam",
        index = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bai"
    output:
        genomecov = temp("results/03_Coverage/{sample}_{aligner}_{markdup}-genome-cov.bed")
    log:
        "results/11_Reports/bedtools/{sample}_{aligner}_{markdup}-genome-cov.log"
    shell:
        "bedtools genomecov "    # Bedtools genomecov, compute the coverage of a feature file among a genome
        "-bga "                   # Report depth in BedGraph format, regions with zero coverage are also reported
        "-ibam {input.markdup} "  # The input file is in BAM format, must be sorted by position
        "1> {output.genomecov} "  # BedGraph output
        "2> {log} "               # Log redirection

###############################################################################
rule samtools_index_markdup:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARKDUP.bam] [INDEX.bai]
    message:
        "SamTools indexing marked as duplicate BAM file {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
       CPUS
    input:
        markdup = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bam"
    output:
        index = "results/02_Mapping/{sample}_{aligner}_{markdup}-mark-dup.bai"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_{markdup}-mark-dup-index.log"
    threads: 
        CPUS
    shell:
        "samtools index "      # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {threads} "        # --threads: Number of additional threads to use (default: 1)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.markdup} "     # Markdup bam input
        "{output.index} "      # Markdup index output
        "&> {log}"             # Log redirection

###############################################################################
rule mark_duplicates_spark:
    # Aim: marking duplicate alignments.
    input:
        calmd = "results/02_Mapping/{sample}_{aligner}_sorted_MD.bam",
    output:
        bam = "results/02_Mapping/{sample}_{aligner}_picard-mark-dup.bam",
        metrics="results/02_Mapping/{sample}_{aligner}_picard-markdup_metrics.txt",
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_picard-mark-dup.log",
    params:
        extra="--remove-sequencing-duplicates",  # optional
        java_opts="",  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.19.1="",  # optional
        #spark_extra="", # optional
    resources:
        mem_mb=16000,
    threads: CPUS
    wrapper:
        "v1.19.1/bio/gatk/markduplicatesspark"

###############################################################################
rule samtools_markdup:
    # Aim: marking duplicate alignments
    # Use: samtools markdup -@ [THREADS] -r -s -O BAM [SORTED.bam] [MARKDUP.bam]
    message:
        "SamTools marking duplicate alignments for {wildcards.sample} sample ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
       CPUS
    input:
        calmd = "results/02_Mapping/{sample}_{aligner}_sorted_MD.bam"
    output:
        markdup = "results/02_Mapping/{sample}_{aligner}_samtools-mark-dup.bam"
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_samtools-mark-dup.log"
    threads: 
        CPUS
    shell:
        "samtools markdup "           # Samtools markdup, tools for alignments in the SAM format with command mark duplicates
        "--threads {threads} "        # -@: Number of additional threads to use (default: 1)
        "-r "                         # -r: Remove duplicate reads
        "-s "                         # -s: Report stats
        "--output-fmt BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.calmd} "              # Sorted bam input
        "{output.markdup} "           # Markdup bam output
        "&> {log}"                    # Log redirection

###############################################################################
rule samtools_calmd:
    # Aim: Generate the MD and the NM tag. The MD tag is for SNP/indel calling without looking at the reference. It does this by carrying information about 
    #the reference that the read does not carry, for a particular alignment. A SNP’s alternate base is carried in the read, but without the MD tag 
    #or use of the alignment reference, it’s impossible to know what the reference base was. Thus, this information is carried in the MD tag. 
    # Use: samtools calmd -@ [THREADS] -b input.bam [REF.FA] > output.bam
    message:
        "SamTools calmd {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
       CPUS,
    params:
        refpath = REFPATH,
        reference = REFERENCE
    input:
        sorted = "results/02_Mapping/{sample}_{aligner}_sorted.bam"
    output:
        calmd = temp("results/02_Mapping/{sample}_{aligner}_sorted_MD.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted_calmd.log"
    threads: 
        CPUS
    shell:
        "samtools calmd "                                               # Samtools calmd, tools to generate the MD tag for SNP/indel calling w/o lookin at the reference
        "--threads {threads} "                                          # -@: Number of additional threads to use (default: 1)
        "-b "                                                           # Output compressed BAM
        "{input.sorted} "                                               # bam input (sorted)
        "{params.refpath}{params.reference} "                           # Reference index filename prefix        
        "1> {output.calmd} "                                            # Sorted bam output
        "2> {log}"                                                      # Log redirection

###############################################################################
rule samtools_sorting:
    # Aim: sorting
    # Use: samtools sort -@ [THREADS] -m [MEM] -T [TMPDIR] -O BAM -o [SORTED.bam] [FIXMATE.bam]
    message:
        "SamTools sorting {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
        CPUS
    resources:
       mem_gb = MEM_GB
    params:
        tmpdir = TMPDIR
    input:
        fixmate = "results/02_Mapping/{sample}_{aligner}_fix-mate.bam"
    output:
        sorted = temp("results/02_Mapping/{sample}_{aligner}_sorted.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sorted.log"
    threads: 
        CPUS
    shell:
        "samtools sort "               # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {threads} "         # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "      # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-T {params.tmpdir} "          # -T: Write temporary files to PREFIX.nnnn.bam
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sorted} "          # Sorted bam output
        "{input.fixmate} "             # Fixmate bam input
        "&> {log}"                     # Log redirection

###############################################################################
rule samtools_fixmate:
    # Aim: filling in mate coordinates
    # Use: samtools fixmate -@ [THREADS] -m -O BAM [SORTBYNAMES.bam] [FIXMATE.bam]
    message:
        "SamTools filling in mate coordinates {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
        CPUS
    input:
        sortbynames = "results/02_Mapping/{sample}_{aligner}_sort-by-names.bam"
    output:
        fixmate = temp("results/02_Mapping/{sample}_{aligner}_fix-mate.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_fixmate.log"
    threads: 
        CPUS
    shell:
        "samtools fixmate "            # Samtools fixmate, tools for alignments in the SAM format with command to fix mate information
        "--threads {threads} "         # -@: Number of additional threads to use (default: 1)
        "-m "                          # -m: Add mate score tag
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.sortbynames} "         # Sortbynames bam input
        "{output.fixmate} "            # Fixmate bam output
        "&> {log}"                     # Log redirection

###############################################################################
rule samtools_sortbynames:
    # Aim: sorting by names
    # Use: samtools sort -t [THREADS] -n -O BAM -o [SORTBYNAMES.bam] [MAPPED.sam]
    message:
        "SamTools sorting by names {wildcards.sample} sample reads ({wildcards.aligner})"
    conda:
        SAMTOOLS
    threads:
        CPUS
    resources:
        mem_gb = MEM_GB
    input:
        mapped = "results/02_Mapping/{sample}_{aligner}-mapped.sam"
    output:
        sortbynames = temp("results/02_Mapping/{sample}_{aligner}_sort-by-names.bam")
    log:
        "results/11_Reports/samtools/{sample}_{aligner}_sort-by-names.log"
    threads: 
        CPUS
    shell:
        "samtools sort "               # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {threads} "         # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "      # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-n "                          # -n: Sort by read name (not compatible with samtools index command)
        "--output-fmt BAM "            # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sortbynames} "     # -o: Write final output to FILE rather than standard output
        "{input.mapped} "              # Mapped reads input
        "&> {log}"                     # Log redirection


###############################################################################
rule bwa_mapping:
    # Aim: reads mapping against reference sequence
    # Use: bwa mem -t [THREADS] -x [REFERENCE] [FWD_R1.fq] [REV_R2.fq] 1> [MAPPED.sam]
    message:
        "BWA-MEM mapping {wildcards.sample} sample reads against reference genome sequence"
    conda:
        BWA
    threads:
        CPUS
    params:
        bwapath = BWAPATH,
        reference = REFERENCE,
        extra = r"-R '@RG\tID:{sample}\tSM:{sample}'"
    input:
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz"
    output:
        mapped = temp("results/02_Mapping/{sample}_bwa-mapped.sam")
    benchmark:
        "benchmarks/bwa/{sample}_bwa-mapped.tsv"
    log:
        "results/11_Reports/bwa/{sample}.log"
    shell:
        "bwa mem "                                                  # BWA-MEM algorithm, performs local alignment.
        "-M "                                                       # Mark shorter split hits as secondary (for Picard compatibility). 
        "-T 0 "                                                     # Don’t output alignment with score lower than INT. This option only affects output.
        "-t {threads} "                                             # -t: Number of threads (default: 12)
        "-v 1 "                                                     # -v: Verbosity level: 1=error, 2=warning, 3=message, 4+=debugging
        "{params.extra} "                                           #
        "{params.bwapath}{params.reference} "                       # Reference index filename prefix
        "{input.fwdreads} "                                         # Forward input reads
        "{input.revreads} "                                         # Reverse input reads
        "1> {output.mapped} "                                       # SAM output
        "2> {log}"                                                  # Log redirection

###############################################################################
rule bowtie2_mapping:
    # Aim: reads mapping against reference sequence
    # Use: bowtie2 -p [THREADS] -x [REFERENCE] -1 [FWD_R1.fq] -2 [REV_R2.fq] -S [MAPPED.sam]
    message:
        "Bowtie2 mapping {wildcards.sample} sample reads against reference genome sequence"
    conda:
        BOWTIE2
    threads:
        CPUS
    params:
        bt2path = BT2PATH,
        reference = REFERENCE,
        sensitivity = SENSITIVITY
    input:
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz"
    output:
        mapped = temp("results/02_Mapping/{sample}_bowtie2-mapped.sam")
    log:
        "results/11_Reports/bowtie2/{sample}.log"
    benchmark:
        "benchmarks/bowtie2/{sample}_bowtie2-mapped.tsv"
    threads: 
        CPUS
    shell:
        "bowtie2 "                                # Bowtie2, an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.
        "--threads {threads} "                    # -p: Number of alignment threads to launch (default: 1)
        "--reorder "                              # Keep the original read order (if multi-processor option -p is used)
        "-x {params.bt2path}{params.reference} "  # -x: Reference index filename prefix (minus trailing .X.bt2) [Bowtie-1 indexes are not compatible]
        "{params.sensitivity} "                   # Preset (default: "--sensitive", same as [-D 15 -R 2 -N 0 -L 22 -i S,1,1.15])
        "-q "                                     # -q: Query input files are FASTQ .fq/.fastq (default)
        "-1 {input.fwdreads} "                    # Forward input reads
        "-2 {input.revreads} "                    # Reverse input reads
        "1> {output.mapped} "                     # -S: File for SAM output (default: stdout)
        "2> {log}"                                # Log redirection

###############################################################################
rule sickle_trim_quality:
    # Aim: windowed adaptive trimming tool for FASTQ files using quality
    # Use: sickle [COMMAND] [OPTIONS]
    message:
        "Sickle-trim low quality sequences trimming for {wildcards.sample} sample"
    conda:
        SICKLETRIM
    params:
        command = COMMAND,
        encoding = ENCODING,
        quality = QUALITY,
        length = LENGTH
    input:
        fwdreads = "results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R1.fastq.gz",
        revreads = "results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R2.fastq.gz"
    output:
        fwdreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R1.fastq.gz",
        revreads = "results/01_Trimming/sickle/{sample}_sickle-trimmed_R2.fastq.gz",
        single = temp("results/01_Trimming/sickle/{sample}_sickle-trimmed_SE.fastq.gz")
    benchmark:
        "benchmarks/sickle/{sample}_sickle-trimmed-qual.tsv"
    log:
        "results/11_Reports/sickle-trim/{sample}.log"
    shell:
       "sickle "                 # Sickle, a windowed adaptive trimming tool for FASTQ files using quality
        "{params.command} "      # Paired-end or single-end sequence trimming
        "-t {params.encoding} "  # --qual-type: Type of quality values, solexa ; illumina ; sanger ; CASAVA, < 1.3 ; 1.3 to 1.7 ; >= 1.8
        "-q {params.quality} "   # --qual-threshold: Threshold for trimming based on average quality in a window (default: 20)
        "-l {params.length} "    # --length-threshold: Threshold to keep a read based on length after trimming (default: 20)
        "-f {input.fwdreads} "   # --pe-file1: Input paired-end forward fastq file
        "-r {input.revreads} "   # --pe-file2: Input paired-end reverse fastq file
        "-g "                    # --gzip-output: Output gzipped files
        "-o {output.fwdreads} "  # --output-pe1: Output trimmed forward fastq file
        "-p {output.revreads} "  # --output-pe2: Output trimmed reverse fastq file (must use -s option)
        "-s {output.single} "    # --output-single: Output trimmed singles fastq file
        "&> {log}"               # Log redirection

###############################################################################
rule cutadapt_adapters_removing:
    # Aim: finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads
    # Use: cutadapt [OPTIONS] -a/-A [ADAPTER] -o [OUT-FWD.fastq.gz] -p [OUT-REV.fastq.gz] [IN-FWD.fastq.gz] [IN-REV.fastq.gz]
    # Rmq: multiple adapter sequences can be given using further -a options, but only the best-matching adapter will be removed
    message:
        "Cutadapt adapters removing for {wildcards.sample} sample"
    conda:
        CUTADAPT
    threads:
        CPUS
    params:
        length = LENGTHc,
        truseq = TRUSEQ,
        nextera = NEXTERA,
        small = SMALL
    input:
        fwdreads = "resources/reads/{sample}_R1.fastq.gz",
        revreads = "resources/reads/{sample}_R2.fastq.gz"
    output:
        fwdreads = temp("results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R1.fastq.gz"),
        revreads = temp("results/01_Trimming/cutadapt/{sample}_cutadapt-removed_R2.fastq.gz")
    log:
        "results/11_Reports/cutadapt/{sample}.log"
    threads: 
        CPUS
    benchmark:
        "benchmarks/cutadapt/{sample}_cutadapt.tsv"
    shell:
       "cutadapt "                           # Cutadapt, finds and removes unwanted sequence from your HT-seq reads
        "--cores {threads} "                 # -j: Number of CPU cores to use. Use 0 to auto-detect (default: 1)
        "--trim-n "                          # --trim-n: Trim N's on ends of reads
        "--minimum-length {params.length} "  # -m: Discard reads shorter than length
        "--adapter {params.truseq} "         # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.truseq} "                # -A: 3' adapter to be removed from second read in a pair
        "--adapter {params.nextera} "        # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.nextera} "               # -A: 3' adapter to be removed from second read in a pair
        "--adapter {params.small} "          # -a: Sequence of an adapter ligated to the 3' end of the first read
        "-A {params.small} "                 # -A: 3' adapter to be removed from second read in a pair
        "--output {output.fwdreads} "        # -o: Write trimmed reads to FILE
        "--paired-output {output.revreads} " # -p: Write second read in a pair to FILE
        "{input.fwdreads} "                  # Input forward reads R1.fastq
        "{input.revreads} "                  # Input reverse reads R2.fastq
        "&> {log}"                           # Log redirection

###############################################################################
rule multiqc_reports_aggregation:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [MULTIQC/]
    priority: 42
    message:
        "MultiQC reports aggregating"
    conda:
        MULTIQC
    input:
        fastqc = "results/00_Quality_Control/fastqc/",
        fastqscreen = "results/00_Quality_Control/fastq-screen/"
    output:
        multiqc = directory("results/00_Quality_Control/multiqc/")
        #report("results/00_Quality/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "results/11_Reports/quality/multiqc.log"
    shell:
        "multiqc "                  # Multiqc, searches in given directories for analysis & compiles a HTML report
        "--quiet "                   # -q: Only show log warning
        "--outdir {output.multiqc} " # -o: Create report in the specified output directory
        "{input.fastqc} "            # Input FastQC files
        "{input.fastqscreen} "       # Input Fastq-Screen
        "--no-ansi "                 # Disable coloured log
        "&> {log}"                   # Log redirection

###############################################################################
rule fastqscreen_contamination_checking:
    # Aim: screen if the composition of the library matches with what you expect
    # Use: fastq_screen [OPTIONS] --outdir [DIR/] [SAMPLE_1.fastq] ... [SAMPLE_n.fastq]
    message:
        "Fastq-Screen reads contamination checking"
    conda:
        FASTQSCREEN
    threads:
        CPUS
    params:
        config = CONFIG,
        mapper = MAPPER,
        subset = SUBSET
    input:
        fastq = "resources/reads/"
    output:
        fastqscreen = directory("results/00_Quality_Control/fastq-screen/")
    log:
        "results/11_Reports/quality/fastq-screen.log"
    shell:
        "fastq_screen "                  # FastqScreen, what did you expect ?
        "-q "                            # --quiet: Only show log warning
        "--threads {threads} "           # --threads: Specifies across how many threads bowtie will be allowed to run
        "--conf {params.config} "        # path to configuration file
        "--aligner {params.mapper} "     # -a: choose aligner 'bowtie', 'bowtie2', 'bwa'
        "--subset {params.subset} "      # Don't use the whole sequence file, but create a subset of specified size
        "--outdir {output.fastqscreen} " # Output director
        "{input.fastq}/*.fastq.gz "      # Input file.fastq
        "&> {log}"                       # Log redirection

###############################################################################
rule fastqc_quality_control:
    # Aim: reads sequence files and produces a quality control report
    # Use: fastqc [OPTIONS] --output [DIR/] [SAMPLE_1.fastq] ... [SAMPLE_n.fastq]
    message:
        "FastQC reads quality controling"
    conda:
        FASTQC
    threads:
        CPUS
    input:
        fastq = "resources/reads/"
    output:
        fastqc = directory("results/00_Quality_Control/fastqc/")
    log:
        "results/11_Reports/quality/fastqc.log"
    shell:
        "mkdir -p {output.fastqc} "   # (*) this directory must exist as the program will not create it
        "2> /dev/null && "            # in silence and then...
        "fastqc "                     # FastQC, a high throughput sequence QC analysis tool
        "--quiet "                    # -q: Supress all progress messages on stdout and only report errors
        "--threads {threads} "        # -t: Specifies files number which can be processed simultaneously
        "--outdir {output.fastqc} "   # -o: Create all output files in the specified output directory (*)
        "{input.fastq}/*.fastq.gz "   # Input file.fastq
        "&> {log}"                    # Log redirection

###############################################################################
