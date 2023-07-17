# SHAVE2: SHort-read Alignment pipeline for VEctor v.2  - with HaplotypeCaller variant caller - SLURM version #

![Developper](<https://badgen.net/badge/Developper/Loïc Talignani/red?scale=0.9>)
![Maintener](<https://badgen.net/badge/Maintener/Loïc Talignani/blue?scale=0.9>)
![MacOSX](<https://badgen.net/badge/icon/Hight Sierra (10.13.6) | Catalina (10.15.7) | Big Sure (11.6.3) | Monterey (12.2.0)/E6055C?icon=apple&label&list=|&scale=0.9>)
![GNU/Linux](<https://badgen.net/badge/icon/Bionic Beaver (18.04) | Focal Fossa (20.04) | Jammy Jellyfish (22.04)/772953?icon=https://www.svgrepo.com/show/25424/ubuntu-logo.svg&label&list=|&scale=0.9>)
![Maintened](<https://badgen.net/badge/Maintened/Yes/red?scale=0.9>)
![Wiki](<https://badgen.net/badge/icon/Wiki/pink?icon=wiki&label&scale=0.9>)
![Open Source](<https://badgen.net/badge/icon/Open Source/purple?icon=https://upload.wikimedia.org/wikipedia/commons/4/44/Corazón.svg&label&scale=0.9>)
![GNU AGPL v3](<https://badgen.net/badge/Licence/GNU AGPL v3/grey?scale=0.9>)
![Github](<https://badgen.net/badge/icon/Github/green?icon=gitlab&label&scale=0.9>)
![Bash](<https://badgen.net/badge/icon/Bash 5.2/black?icon=terminal&label&scale=0.9>)
![Python](<https://badgen.net/badge/icon/Python 3.10.6/black?icon=https://upload.wikimedia.org/wikipedia/commons/0/0a/Python.svg&label&scale=0.9>)
![Snakemake](<https://badgen.net/badge/icon/Snakemake 6.12.1/black?icon=https://upload.wikimedia.org/wikipedia/commons/d/d3/Python_icon_%28black_and_white%29.svg&label&scale=0.9>)
![Conda](<https://badgen.net/badge/icon/Conda 4.10.3/black?icon=codacy&label&scale=0.9>)


## ~ ABOUT ~ ##

SHAVE2 is a bioinformatic pipeline used for mosquitoes (*Aedes / Anopheles*) genome alignments from Illumina short reads, based on GATK Best Practices, except for the BQSR and VQSR steps. 

In brief, SHAVE2 remove adatpers, report quality reads, aligns reads to a reference genome, fix incorrect mates, mark duplicates, add indel qualities to BAM files and call variants and genotypes.

### Note about Variant and Genotype Calling ###

Indel realignment was dropped by Broad Institute about three years ago, as they found that this step was no longer useful when the variant calling was done with HaplotypeCaller or Mutect2, which implement a more sophisticated and effective form of realignment. [see here](https://gatk.broadinstitute.org/hc/en-us/articles/360038154292/comments/360010778692). SHAVE2 use HaplotypeCaller in ERC mode to call variants and obtain the genotype likelihoods. Then, GenomicsDBImport import single-sample GVCFs into GenomicsDB before GenotypeGVCFs perform a joint genotyping on samples pre-called with HaplotypeCaller. The last major step is a fully customizable hard-filtering, using GATK VariantFiltration instead of Variant Qualtity Score Recalibration (VQSR).

### Note about BQSR and VQSR ###

*Base Quality Score Recalibration* step needs as input a *known variation* VCF file, refering to the Ensembl-Variation database or dbSNP database who stores areas of genome that differ between individual genomes ("variants"). However, we do not have any prior list of known variants for our mosquito species, that's why we cannot do BQSR. 

But the value of BQSR is [increasingly being questioned](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1279-z#Abs1) as mappers and callers are typically updated. As an example: using HaplotypeCaller instead of UnifiedGenotyper greatly improves the handling of indels.

*Variant Quality Score recalibration* is probably the hardest part of the Best Practices to get right [according to Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/4413056319131-VariantRecalibrator). In a nutshell, it is a sophisticated filtering technique applied on the variant callset that uses machine learning to model the technical profile of variants in a training set and uses that to filter out probable artifacts from the callset.

The key point is that it use known, highly validated variant resources (omni, 1000 Genomes, hapmap) to select a subset of variants within our callset that we’re really confident are probably true positives (the training set). Unfortunately, no highly validated variant resource is available for mosquitoes at this time, so we decided to apply hard-filtering and leave the choice of parameters to the user. 

Written for **MOVE-ADAPT** project. 

### Features ###
- Control reads quality (_multiQC html report_) and clean it  
- Align reads (_bam files_), 
- Mark duplicates to BAM files,
- Add MD and NM tag to BAM files,
- Fix mates,
- BAM file validation according to SAM/BAM specifications,
- Variants calling (_vcf files_),
- Genotyping, 
- VCF compression,
- VCF indexing

### Version ###
*V2.2023.07.13*  

### Directed Acyclic Graph ###
<img src="./visuals/dag.png" width="200" height="400">


## ~ INSTALLATIONS ~ ##


# SHAVE2 #
Clone _(HTTPS)_ the [SHAVE2](https://github.com/ltalignani/SHAVE2) repository on github:

```shell
git clone git@github.com:ltalignani/SHAVE2.git

```

Difference between **Download** and **Clone**:  
- To create a copy of a remote repository’s files on your computer, you can either **Download** or **Clone** the repository  
- If you download it, you **cannot sync** the repository with the remote repository on GitHub  
- Cloning a repository is the same as downloading, except it preserves the Git connection with the remote repository  
- You can then modify the files locally and upload the changes to the remote repository on GitHub  
- You can then **update** the files locally and download the changes from the remote repository on GitHub  
```shell
git pull --verbose
```
Then, upload the SHAVE2/ archive on your slurm cluster. 

## ~ USAGE ~ ##

1. Copy your **paired-end** reads in **.fastq.gz** format files into: **.raw/** directory
2. Execute the **sbatch Start_shave2_slurm.sh** bash script to run the SHAVE2 pipeline

```shell
sbatch Start_shave2_slurm.sh
```
Yours analyzes will start with default configuration settings  

_Option-1: Edit **config.yaml** file in **./config/** directory_  
_Option-2: Edit **fastq-screen.conf** file in **./config/** directory_  


## ~ RESULTS ~ ##

Yours results are available in **./results/** directory, as follow:  
_(file names keep track which tools / params was used: \<**sample**\>\_\<**aligner**\>\_\<**mincov**\>)_  

### 00_Quality_Control ###
| File | Object |
|:--- | :--- |
| **fastq-screen** | raw reads putative contaminations reports for each samples, in _html_, _png_ and _txt_ formats |
| **fastqc** | raw reads quality reports for each samples, in _html_ and _zip_ formats |
| **multiqc** | fastq-screen and fastqc results agregation report for all samples, in _html_ format |
| **qualimap** | facilitate the quality control of alignment sequencing data and its derivatives like feature counts, in _html_ format |
| **validatesamfile** | This tool reports on the validity of a SAM or BAM file relative to the SAM format specification, in _txt_ format |

### 01_Trimming ###
| File | Object |
|:--- | :--- |
| **trimmomatic/ directory** | paired reads, without adapters and quality trimmed, in _fastq.gz_ format

### 02_Mapping ###
| File | Object |
|:--- | :--- |
| **mark-dup.bam** | read alignments, in _bam_ format _(can be visualized in, i.e. IGV)_ |
| **mark-dup.bai** | bam indexes _bai_ use in i.e. IGV with _./resources/genomes/AalbF3.fasta_ |
| **markdup_metrics.txt** | bam metrics created by Picard MarkDuplicates |

### 04_Polishing ###
| File | Object |
|:--- | :--- |
| **md_fixed.bam** | polished bam files , in _bam_ format |
| **md_fixed.bam.bai** | polished bam indexes, in _bai_ format |

### 05_Variants ###

| File | Object |
|:--- | :--- |
| _variant-call.g.vcf_ | SNVs and Indels calling in _gvcf_ format | 
| _variant-call.g.vcf.idx_ | SNV and Indels calling in _idx_ format |  
| _genomicdb_ | g.vcfs merged by GenomicsDBImport for GenotypeGVCFs| 
| _genotyped.vcf_ | SNVs genotyped in _vcf_ format |. 
| _genotyped.vcf.idx_ | genotype index in _idx_ format (automatically generated by GenotypeGVCFs) |. 
| _combinedGVCF.vcf.gz_ | alternative to GenomicsDBImport, is a merging of GVCFs that can eventually be input into GenoTypeGVCFs |
| _GenotypeGVCFs.hf.vcf.gz_ | genotyped SNVs hard-filtered, in _vcf_ format (default config: tempdir, removed, save disk usage)_ |

### 10_graphs ###
| File | Object |
|:--- | :--- |
| **dag** | directed acyclic graph of jobs, in _pdf_ and _png_ formats |
| **rulegraph** | dependency graph of rules, in _pdf_ and _png_ formats _(less crowded than above DAG of jobs, but also show less information)_  |
| **filegraph** | dependency graph of rules with their input and output files in the dot language, in _pdf_ and _png_ formats _(an intermediate solution between above DAG of jobs and the rule graph)_ |

### 11_Reports ###
- All _non-empty_ **log** for each tool and each sample
- files_summary.txt: summary of all files created by the workflow, in _txt_ format  
_(columns: filename, modification time, rule version, status, plan)_


## ~ CONFIGURATION ~ ##

If needed, see or edit default settings in **config.yaml** file in **./config/** directory  

### RESOURCES ###
**samples** and **units** aren't necessary for the moment (update in progress)

### MODULES ###
Modules that will be loaded for each tool. Change the path to correspond to your cluster. 

### REFERENCE ###
- **ref_name**: reference name used for genome mapping _(default config: 'Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4')_
- **path**: path to genome reference
- **reference**: reference sequence path used for genome mapping _(default config: 'AGAMP4')_
- **index**: reference sequence index used for genome mapping _(default config: 'Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.fai')_
- **dictionary**: reference sequence dictionary used by GATK's tools _(default config: 'Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.dict')_

### PROCESSING ###
- **tmpdir**: temporary directory _(default config: '$TMPDIR')_

### GATK ###
- **genomicsdbimport**: genomicsdbimport flags _(default config: "")_
- **genotypegvcfs**: genotypegvcfs flags _(default config: '-include-non-variant-sites')_
- **haplotypecaller**: haplotypecaller flags _(default config: '-ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES')_

### KNOWN VARIATION SITES ###
- **alleles_target**: path to VCF containing known varition sites _(default config: "--alleles xxxx")_. Replace "xxxx" by path to the VCF containing known variation site.

### BWA###
- **path**: path the BWA indexes _(default config: 'resources/indexes/bwa/')_

### FASTQ-SCREEN ###
- **config**: path to the fastq-screen configuration file _(default config: ./config/fastq-screen.conf)_
- **subset**: do not use the whole sequence file, but create a temporary dataset of this specified number of read _(default config: '1000')_
- **aligner**: specify the aligner to use for the mapping. Valid arguments are 'bowtie', bowtie2' or 'bwa' _(default config: 'bwa')_

### VARANTFILTRATION ###
- **vqsr**: not coded at this time
- **hard**: settings for VCFs hardfiltering _(default config:myfilter": "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0)_


### Glossary ###
- **BAM**: Binary Alignment Map
- **BAI**: BAM Indexes
- **FASTA**: Fast-All
- **FASTQ**: FASTA with Quality
- **FAI**: FASTA Indexes
- **SAM**: Sequence Alignment Map

### Directories tree structure ###
```shell
🖥️️ Start_shave2_slurm.sh
📚 README.md
📚 LICENSE
📂 Cluster_logs/
📂 visuals/
 └── 📈 DAG.png
📂 config/
 ├── ⚙️ config.yaml
 └── ⚙️ fastq-screen.conf
📂 raw/
 │    ├── 🛡️ .gitkeep
 │    ├── 📦 ERR3343471_R1.fastq.gz
 │    └── 📦 ERR3343471_R2.fastq.gz
📂 resources/
 ├── 📂 genomes/
 │    └── 🧬 AalbF3.fasta
 ├── 📂 indexes/
 │    ├── 📂 bwa/
 │    │    └── 🗂️ AGAMP4
 ├── 📂 Adapters
 │    ├── 🧬 nextera
 │    ├── 🧬 truseq2-pe
 │    ├── 🧬 truseq2-se
 │    ├── 🧬 truseq3-pe
 │    ├── 🧬 truseq3-pe2
 │    └── 🧬 truseq3-se  
 │
 ├── 📂 results/
 ├── 📂 slurm/
 │    ├── ⚙️ config.yaml
 │    └── 🖥️️ status-sacct.sh   
📂 workflow/
 ├── 📂 envs/
 │    ├── 🍜 bcftools-1.14.yaml
 │    ├── 🍜 bedtools-2.30.0.yaml
 │    ├── 🍜 bowtie2-2.4.4.yaml
 │    ├── 🍜 bwa-0.7.17.yaml
 │    ├── 🍜 cutadapt-3.5.yaml
 │    ├── 🍜 fastq-screen-0.14.0.yam
 │    ├── 🍜 fastqc-0.11.9.yaml
 │    ├── 🍜 gatk-3.8.yaml
 │    ├── 🍜 gatk-4.3.0.0.yaml
 │    ├── 🍜 gawk-5.1.0.yaml
 │    ├── 🍜 lofreq-2.1.5.yaml
 │    ├── 🍜 multiqc-1.11.yaml
 │    ├── 🍜 picard-2.27.4.yaml
 │    ├── 🍜 samtools-1.14.yaml
 │    └── 🍜 sickle-trim-1.33.yaml
 ├── 📂 report/
 ├── 📂 rules/
 │    ├── 📜 calling.smk
 │    ├── 📜 common.smk
 │    ├── 📜 filtering.smk
 │    ├── 📜 mapping.smk
 │    ├── 📜 polishing.smk
 │    ├── 📜 stats.smk
 │    └── 📜 vcf_stats.smk
 ├── 📂 schemas/     
 ├── 📂 scripts/
 │    ├── 🍜 report_vcf.rmd
 │    └── 🍜 report.rmd
 └── 📜 snakefile.smk
```

## ~ SUPPORT ~ ##
1. Read The Fabulous Manual!
2. Read de Awsome Wiki! (todo...)
3. Create a new issue: Issues > New issue > Describe your issue
4. Send an email to [loic.talignani@ird.fr](url)

## ~ ROADMAP ~ ##
- Add a wiki!  

## ~ AUTHORS & ACKNOWLEDGMENTS ~ ##
- Loïc TALIGNANI (Developer and Maintener)  
- I would like to thanks **Sebastien Ravel**, CIRAD ingeneer, for his great help for the variant calling part.


## ~ CONTRIBUTING ~ ##
Open to contributions!  
Testing code, finding issues, asking for update, proposing new features...  
Use Git tools to share!  

## ~ PROJECT STATUS ~ ##
This project is **regularly updated** and **actively maintened**  
However, you can be volunteer to step in as **developer** or **maintainer**  

For information about main git roles:  
- **Guests** are _not active contributors_ in private projects, they can only see, and leave comments and issues
- **Reporters** are _read-only contributors_, they can't write to the repository, but can on issues  
- **Developers** are _direct contributors_, they have access to everything to go from idea to production  
_Unless something has been explicitly restricted_  
- **Maintainers** are _super-developers_, they are able to push to master, deploy to production  
_This role is often held by maintainers and engineering managers_  
- **Owners** are essentially _group-admins_, they can give access to groups and have destructive capabilities  

## ~ LICENSE ~ ##
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)  

## ~ REFERENCES ~ ##
**Sustainable data analysis with Snakemake**  
Felix Mölder, Kim Philipp Jablonski, Brice Letcher, Michael B. Hall, Christopher H. Tomkins-Tinch, Vanessa Sochat, Jan Forster, Soohyun Lee, Sven O. Twardziok, Alexander Kanitz, Andreas Wilm, Manuel Holtgrewe, Sven Rahmann, Sven Nahnsen, Johannes Köster  
_F1000Research (2021)_  
**DOI**: [https://doi.org/10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)  
**Publication**: [https://f1000research.com/articles/10-33/v1](https://f1000research.com/articles/10-33/v1)  
**Source code**: [https://github.com/snakemake/snakemake](https://github.com/snakemake/snakemake)  
**Documentation**: [https://snakemake.readthedocs.io/en/stable/index.html](https://snakemake.readthedocs.io/en/stable/index.html)  

**Tabix: fast retrieval of sequence features from generic TAB-delimited files**  
Heng Li  
_Bioinformatics, Volume 27, Issue 5 (2011)_  
**DOI**: [https://doi.org/10.1093/bioinformatics/btq671](https://doi.org/10.1093/bioinformatics/btq671)  
**Publication**: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/)  
**Source code**: [https://github.com/samtools/samtools](https://github.com/samtools/samtools)  
**Documentation**: [http://samtools.sourceforge.net](http://samtools.sourceforge.net)  

**GATK: A MapReduce framework for analyzing next-generation DNA sequencing data**
_Genome Research, Volume 20: 1297-1303 (2010)_
**DOI**: [https://doi.org/10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110)
**Publication**: [https://genome.cshlp.org/content/20/9/1297](https://genome.cshlp.org/content/20/9/1297)
**Source code**:[https://github.com/broadinstitute/gatk](https://github.com/broadinstitute/gatk)
**Documentation**:[https://gatk.broadinstitute.org/hc/en-us](https://gatk.broadinstitute.org/hc/en-us)

**Picard-tools: **
_Broad Institute, GitHub repository (2019)_
**DOI**:
**Publication**: 
**Source code**:https://github.com/broadinstitute/picard](https://github.com/broadinstitute/picard)
**Documentation**:[https://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)

**The AWK Programming Language**  
Al Aho, Brian Kernighan and Peter Weinberger  
_Addison-Wesley (1988)_  
**ISBN**: [https://www.biblio.com/9780201079814](https://www.biblio.com/9780201079814)  
**Publication**: []()  
**Source code**: [https://github.com/onetrueawk/awk](https://github.com/onetrueawk/awk)  
**Documentation**: [https://www.gnu.org/software/gawk/manual/gawk.html](https://www.gnu.org/software/gawk/manual/gawk.html)  

**Twelve years of SAMtools and BCFtools**  
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies and Heng Li  
_GigaScience, Volume 10, Issue 2 (2021)_  
**DOI**: [https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)  
**Publication**: [https://academic.oup.com/gigascience/article/10/2/giab008/6137722](https://academic.oup.com/gigascience/article/10/2/giab008/6137722)  
**Source code**: [https://github.com/samtools/samtools](https://github.com/samtools/samtools)  
**Documentation**: [http://samtools.sourceforge.net](http://samtools.sourceforge.net)  

**Fast and accurate short read alignment with Burrows-Wheeler Transform**  
Heng Li and Richard Durbin  
_Bioinformatics, Volume 25, Aricle 1754-60 (2009)_  
**DOI**: [https://doi.org/10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324)  
**Publication**: [https://pubmed.ncbi.nlm.nih.gov/19451168@](https://pubmed.ncbi.nlm.nih.gov/19451168)  
**Source code**: [https://github.com/lh3/bwa](https://github.com/lh3/bwa)  
**Documentation**: [http://bio-bwa.sourceforge.net](http://bio-bwa.sourceforge.net)  

**Trimmomatic: A sliding-window, adaptive, quality-based trimming tool for FastQ files**  
Joshi NA and Fass JN  
_(2011)  
**DOI**: [https://doi.org/10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)  
**Publication**: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/pdf/btu170.pdf](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/pdf/btu170.pdf)  
**Source code**: [https://github.com/usadellab/Trimmomatic](https://github.com/usadellab/Trimmomatic)  
**Documentation**: [https://github.com/usadellab/Trimmomatic](https://github.com/usadellab/Trimmomatic)  

**MultiQC: summarize analysis results for multiple tools and samples in a single report**  
Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller  
_Bioinformatics, Volume 32, Issue 19 (2016)_  
**DOI**: [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)  
**Publication**: [https://academic.oup.com/bioinformatics/article/32/19/3047/2196507](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)  
**Source code**: [https://github.com/ewels/MultiQC](https://github.com/ewels/MultiQC)  
**Documentation**: [https://multiqc.info](https://multiqc.info)  

**FastQ Screen: A tool for multi-genome mapping and quality control**  
Wingett SW and Andrews S  
_F1000Research (2018)_  
**DOI**: [https://doi.org/10.12688/f1000research.15931.2](https://doi.org/10.12688/f1000research.15931.2)  
**Publication**: [https://f1000research.com/articles/7-1338/v2](https://f1000research.com/articles/7-1338/v2)  
**Source code**: [https://github.com/StevenWingett/FastQ-Screen](https://github.com/StevenWingett/FastQ-Screen)  
**Documentation**: [https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen)  

**FastQC: A quality control tool for high throughput sequence data**  
Simon Andrews  
_Online (2010)_  
**DOI**: [https://doi.org/](https://doi.org/)  
**Publication**: []()  
**Source code**: [https://github.com/s-andrews/FastQC](https://github.com/s-andrews/FastQC)  
**Documentation**: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)  
