#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/virtual_spike_in/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/virtual_spike_in/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu
# gen_featurecounts.sbatch --- Generate Featurecounts
#
# Filename: gen_featurecounts.sbatch
# Description: Make Featurcounts with Coordinate Modification
# Author: Student Zachary Maas <zama8258@colorado.edu>
# Maintainer: Student Zachary Maas <zama8258@colorado.edu>
#

# Commentary:
#
# This file contains code for generating counts tables using
# featureCounts, with builtin support for coordinate modification for
# modifying regions in the bedfile.
#

# Code:

# Handle Errors Properly:
# -e - exits on the first error
# -u - fails when unset variables are called
# -o pipefail - fails when error in a pipe
set -euo pipefail

# Recently tweaked for better POSIX compliance.
function logr() {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Dynamically load bedtools if it isn't available
if ! type -t bedtools; then
		module load bedtools
fi

# Set the number of available cores
NUM_CORES=8

# TmpDir=/scratch/Users/zama8258/processed_nascent/scratch/features
CramDir=/scratch/Shares/dowell/ripin/dm6_mapped_PRO-seq/mapped/crams
BamDir=/scratch/Users/zama8258/virtual_spike_in/scratch/counts/bams
BaseDir=/scratch/Users/zama8258/virtual_spike_in/scratch
# Infile="$BaseDir"/fpkm/LacI_plus_Cu1_S1_R1_001.sorted.isoform_max.bed
Infile=/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed

outFull="$BaseDir"/counts/counts_full.txt

safFull="$BaseDir"/counts/full.saf

logr "Adjusting Coordinates"
awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$Infile" > "$safFull"

## Change directory because FeatureCounts is picky about running in
## the same directory as the crams.
logr "Changing Directories"
# pushd "$CramDir" || exit
# logr "Generating Bam Files"

module load samtools
CRAM_REFERENCE=/scratch/Shares/dowell/genomes/dm6/dm6.fa
# find . -maxdepth 1 -name '*.cram' | \
		# 		sed -e 's/.cram//g' | \
		# 		parallel samtools view -b -T "$CRAM_REFERENCE" -o "$BamDir"/{}.bam {}.cram

# popd
pushd "$CramDir" || exit
# parallel samtools index ::: *.bam
logr "Intersecting the Reference Sequence"
/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		-F 'SAF' \
		-a "$safFull" \
		-o "$outFull" \
		./*.bam

tail -n +2 "$outFull" > "$outFull"_without_header

popd

logr "Done"

#
# gen_featurecounts.sbatch ends here
