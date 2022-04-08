#!/bin/bash
#SBATCH --job-name=nextflow # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=zama8258@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=52:00:00 # Time limit hrs:min:sec
#SBATCH --partition long
#SBATCH --mem=8gb # Memory limit

################### SET VARIABLES ######################################

export PATH=~:$PATH
export PATH=~/.local/bin:$PATH
MAIN=/scratch/Users/zama8258/Nascent-Flow/main.nf

########################################################################
################### LOAD NECESSARY MODULES #############################

module load sra/2.8.0
module load samtools/1.8
module load hisat2/2.1.0
module load bedtools/2.25.0
module load gcc/7.1.0
module load seqkit/0.9.0
module load fastqc/0.11.8
module load bbmap/38.05
module load igvtools/2.3.75
module load preseq/2.0.3
module load mpich/3.2.1
module load python/2.7.14/rseqc/2.6.4

########################################################################
################## PRINT JOB INFO ######################################

printf "Sample ID: $ROOTNAME"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################

nextflow run $MAIN -profile dm6 \
         --genome_id dm6 --singleEnd -resume \
         --sras '/scratch/Users/zama8258/virtual_spike_in/dukler_analysis/drosophila/sra/*.sra' \
         --workdir /scratch/Users/zama8258/work_gro/ \
         --email zachary.maas@colorado.edu \
         --outdir /scratch/Users/zama8258/virtual_spike_in/dukler_analysis/drosophila/nascent_flow_out/
