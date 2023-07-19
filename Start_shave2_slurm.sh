#!/bin/bash
###################configuration slurm##############################
#SBATCH -A aedes_amplicon
#SBATCH --job-name=shave2
#SBATCH --time=6-23:00:00
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=1GB
#SBATCH -o Cluster_logs/shave2-%x-%j-%N.out
#SBATCH -e Cluster_logs/shave2-%x-%j-%N.err
#SBATCH --mail-user=loic.talignani@ird.fr
#SBATCH --mail-type=ALL
###################################################################

# USAGE: sbatch Start_unifiedgenotyper.sh

###Charge module
module purge
module load snakemake/7.25.0
module load graphviz/2.40.1
module load r/4.2.1

echo -e "${green} Current directory is $(pwd) ${nc}"

# set umask to avoid locking each other out of directories
umask 002

##### Colors ######
red="\033[1;31m"   # red
green="\033[1;32m" # green
ylo="\033[1;33m"   # yellow
blue="\033[1;34m"  # blue
nc="\033[0m"       # no color

###### About ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}#####${nc} ${red}ABOUT${nc} ${green}#####${nc}"
echo -e "${green}-----------------${nc}"
echo ""
echo -e "${blue}Name${nc} __________________ Start_shave2_slurm.sh"
echo -e "${blue}Author${nc} ________________ Loïc Talignani"
echo -e "${blue}Affiliation${nc} ___________ UMR_MIVEGEC"
echo -e "${blue}Aim${nc} ___________________ Bash script for ${red}SH${nc}ort-read ${red}A${nc}lignment pipeline for ${red}VE${nc}ctor v.2"
echo -e "${blue}Date${nc} __________________ 2023.07.18"
echo -e "${blue}Run${nc} ___________________ bash Start_shave2_slurm.sh"
echo -e "${blue}Latest Modification${nc} ___ Updated for IFB-core"


# Set working directory
workdir=$(pwd)            #$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
max_threads="30"

echo -e "Workdir is "${workdir}

###### Rename samples ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}##############${nc} ${red}RENAME FASTQ FILES${nc} ${green}##############${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

# Rename fastq files to remove "_001" Illumina pattern.
## De/comment (#) if you want keep Illumina barcode-ID and/or Illumina line-ID
#rename "s/_S\d+_/_/" ${workdir}/raw/*.fastq.gz                  # Remove barcode-ID like {_S001_}
rename "s/_L\d+_/_/" ${workdir}/raw/*.fastq.gz                  # Remove line-ID ID like {_L001_}
rename "s/_001.fastq.gz/.fastq.gz/" ${workdir}/raw/*.fastq.gz   # Remove end-name ID like {_001}.fastq.gz
#rename 's/(\w_).*_(R[1-2]).*(.fastq.gz)/$1$2$3/' *.fastq.gz                 # Keep only expr. in ( )

###### Call snakemake pipeline ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE START${nc} ${green}###########${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

echo -e "${blue}Unlocking working directory:${nc}"
echo ""

snakemake --profile slurm/ --directory ${workdir}/ --snakefile workflow/snakefile.smk --unlock

echo ""
echo -e "${blue}Let's run!${nc}"
echo ""

snakemake --profile slurm/ --snakefile workflow/snakefile.smk --cores 30 --configfile config/config.yaml

###### Create usefull graphs, summary and logs ######
echo ""
echo -e "${green}------------------------------------------------------------------------${nc}"
echo -e "${green}###########${nc} ${red}SNAKEMAKE PIPELINE LOGS${nc} ${green}############${nc}"
echo -e "${green}------------------------------------------------------------------------${nc}"
echo ""

mkdir ${workdir}/results/10_Graphs/ 2> /dev/null

graph_list="dag rulegraph filegraph"
extention_list="pdf png"

for graph in ${graph_list} ; do
    for extention in ${extention_list} ; do
	snakemake --profile slurm/ --snakefile workflow/snakefile.smk --${graph} | dot -T${extention} > ${workdir}/results/10_Graphs/${graph}.${extention} ;
    done ;
done

snakemake --directory ${workdir} --profile slurm/ --snakefile workflow/snakefile.smk --summary > ${workdir}/results/11_Reports/files_summary.txt

###### End managment ######
echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo -e "${blue}##################${nc} ${red}SCRIPT END${nc} ${blue}###################${nc}"
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""

find ${workdir}/results/ -type f -empty -delete                 # Remove empty file (like empty log)
find ${workdir}/results/ -type d -empty -delete                 # Remove empty directory

time_stamp_end=$(date +"%Y-%m-%d %H:%M")                        # Get date / hour ending analyzes
elapsed_time=${SECONDS}                                         # Get SECONDS counter
minutes=$((${elapsed_time}/60))                                 # / 60 = minutes
seconds=$((${elapsed_time}%60))                                    # % 60 = seconds

echo -e "${red}End Time${nc} ______________ ${time_stamp_end}"                                                             # Print analyzes ending time
echo -e "${red}Processing Time${nc} _______ ${ylo}${minutes}${nc} minutes and ${ylo}${seconds}${nc} seconds elapsed"       # Print total time elapsed

echo ""
echo -e "${blue}------------------------------------------------------------------------${nc}"
echo ""
