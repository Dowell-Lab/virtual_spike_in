#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/virtual_spike_in/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/virtual_spike_in/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu

set -euxo pipefail

# Scripts
fc=/scratch/Users/zama8258/subread-1.6.2-Linux-x86_64/bin/featureCounts
NUM_CORES=8

# IO Vars
refGenes=/scratch/Users/zama8258/virtual_spike_in/scratch/fpkm/human/dG3BP_1hAs_U2OS_b1.sorted.isoform_max.bed
bamDir=/scratch/Shares/dowell/ripin/PRO-seq/mapped/bams
outDir=/scratch/Users/zama8258/virtual_spike_in/scratch/norm
baseDir="$outDir"
outBed="$outDir"/coords.bed

# Filter to things longer than 120kb  and adjust the genomic coordinates
awk -v OFS='\t' -f <(cat - <<- "EOF"
BEGIN {}
{
if ($3 - $2 > 120501) {
  if ($6 == "+") {
    print $1, $2 + 120000, $3 - 500, $4, $5, $6
  } else {
    print $1, $2 + 500, $3 - 120000, $4, $5, $6
  }
}
}
END {}
EOF
										) "$refGenes" > "$outBed"


# Generate featurecounts
mkdir -p "$baseDir"/counts/
outFile="$baseDir"/counts/counts_long_ends.txt
safFile="$baseDir"/counts/full.saf
Infile="$outBed"
awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$Infile" > "$safFile"

bams="$(find "$bamDir" -name "*.bam" | xargs echo)"
cd "$bamDir" || exit
"$fc" \
		-T "$NUM_CORES" \
		-s 1 \
		-F 'SAF' \
		-a "$safFile" \
		-o "$outFile" \
		$bams
