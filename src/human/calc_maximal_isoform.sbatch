#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/virtual_spike_in/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/virtual_spike_in/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=32gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zama8258@colorado.edu

## This script contains code to filter out the maximal isoform of all
## genes in a bedfile using bedops bedmap.

set -euxo pipefail

## Start timing the script
SECONDS=0

## This is a logrging function that will display the elapsed time of
## the script and how long the script has been executing for on the
## server. Useful for tracking execution time in various scripts.
function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Make sure we have the necessary modules on the cluster
if ! type -t bedtools
then module load bedtools
fi


##############################
## Variables we always need ##
##############################

# This is the file we're interested in analyzing
# InterestFile=$0 ## TODO	argument parsing
#InterestFile=/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/$srr.tri.BedGraph
RefSeq=/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed

## This determines the directory that we are executing the current
## script from, so we can call other scripts in the	same directory.
ScriptDir=/scratch/Users/zama8258/pause_analysis_src

## This is a scratch directory that we use to hold temporary files.
ScratchDir=/scratch/Users/zama8258/virtual_spike_in/scratch/fpkm/human

## This is the output directory
OutDir="$ScratchDir"

#####################################
## Generate	Temporary Files in RAM ##
#####################################

## Here we generate the temporary directory used to store the temp
## files we use during file processing. This is a fast process because
## we're only reading and writing to RAM, which gives us IO on the
## order of 60Gb/s
logr "Running ""$InterestFile"" in Production Mode"
TmpDir=$(mktemp -d)

PosRefFile="$TmpDir""/""$(uuidgen)"
NegRefFile="$TmpDir""/""$(uuidgen)"

RootName=$(basename "$InterestFile" .bam)

# Cache the sum files	because they take the longest to generate
PosStrandSums="$ScratchDir"/"$RootName".pos.count
NegStrandSums="$ScratchDir"/"$RootName".neg.count

PosCountsMapped=0
NegCountsMapped=0

PosScalingFactor=0
NegScalingFactor=0

PosStrandMillionsNormalizedSums="$TmpDir""/""$(uuidgen)"
NegStrandMillionsNormalizedSums="$TmpDir""/""$(uuidgen)"

PosStrandFpkmNormalizedSums="$TmpDir""/""$(uuidgen)"
NegStrandFpkmNormalizedSums="$TmpDir""/""$(uuidgen)"

PosStrandFinalSorted="$ScratchDir"/"$RootName".pos.maxiso
NegStrandFinalSorted="$ScratchDir"/"$RootName".neg.maxiso

StrandsMerged="$ScratchDir"/"$RootName".merge.maxiso
FPKMCommonID="$ScratchDir"/"$RootName".merge.maxiso.common
FinalFPKM="$OutDir"/"$RootName".isoform_max
FinalOut="$FinalFPKM".bed

# Clean up temp files on exit
function cleanup {
		rm -rf "$TmpDir"
		logr "Deleted temporary directory $TmpDir"
}
# Register the cleanup function to be called on the EXIT signal
trap cleanup EXIT

###############################
## Start Performing Analysis ##
###############################

## Here, we perform strand-specific filtering to account for `bedtools
## multicov` being unable to write out separate results for each
## strand. This strips down the RefSeq annotation to a bed6 format,
## and drops out any non NM_ coded genes, since for the purpose of
## this script we only care about protein coding genes.
logr "Filtering RefSeq Annotations by Strand for ""$InterestFile"	&
grep NM "$RefSeq" | awk -v OFS='\t' '{if ($6 == "+") print $1, $2, $3, $4, $5, $6}'\
												> "$PosRefFile" &
grep NM "$RefSeq" | awk -v OFS='\t' '{if ($6 == "-") print $1, $2, $3, $4, $5, $6}'\
												> "$NegRefFile" &
wait

## Here, we make sure that our bam file is indexed and if not, perform
## the indexing
# if [ ! -f "$InterestFile".bai ]; then
# 		logr "No Index File found. Generating One."
# 		module load samtools
# 		samtools index "$InterestFile"
# fi

if [ ! -f "$PosStrandSums" ] || [ ! -f "$NegStrandSums" ]; then
		## Here, we construct a bedfile from our input BAM, in order to
		## normalize it by FPKM. We do this using bedtools multicov on the
		## sample in a strand specific manner. Note that bedtools multicov
		## will not account for duplicate reads by default, so if we want RPKM
		## instead of FPKM we need to add	the	`-D` flag
		export CRAM_REFERENCE=/scratch/Shares/dowell/genomes/dm6/dm6.fa
		logr "Calculating Strand-Specific Coverage for BAM" &
		bedtools multicov	-s -bams "$InterestFile" -bed "$PosRefFile" > "$PosStrandSums" &
		bedtools multicov	-s -bams "$InterestFile" -bed "$NegRefFile" > "$NegStrandSums"	&
		wait
fi
## TODO This is where the error originates...

## Next, we need to find the sum of all reads mapped in our regions of
## interest, for FPKM. We do this by adding up every mapped count in
## the 5th column of our bedfile from bedtools multicov
read -r PosCountsMapped <<< "$(awk	-F '\t' -v OFS='\t'	'{sum += $7} END {print sum}' "$PosStrandSums")"
read -r NegCountsMapped <<< "$(awk	-F '\t' -v OFS='\t'	'{sum += $7} END {print sum}' "$NegStrandSums")"

## Then, we generate a scaling factor for RPKM by dividing our mapped
## counts by 1 million
read -r PosScalingFactor <<< "$(echo "$PosCountsMapped" / 1000000 | bc -l)"
read -r NegScalingFactor <<< "$(echo "$NegCountsMapped" / 1000000 | bc -l)"

echo "$PosScalingFactor" "$NegScalingFactor"

## We proceed to divide every count in our bedfile by the
## strand-specific scaling factor
awk	-F '\t' -v OFS='\t' -v norm="$PosScalingFactor" \
		'{print $1, $2, $3, $4, $5, $6, $7 / norm}' "$PosStrandSums" > "$PosStrandMillionsNormalizedSums"
awk	-F '\t' -v OFS='\t' -v norm="$NegScalingFactor" \
		'{print $1, $2, $3, $4, $5, $6, $7 / norm}' "$NegStrandSums" > "$NegStrandMillionsNormalizedSums"

## We finally take the sums (now normalized for millions mapped), then
## divide those values by gene length.
awk -F '\t' -v OFS='\t'	'{print $1, $2, $3, $4, $5, $6, $7 / ($3 - $2)}' \
		"$PosStrandMillionsNormalizedSums" >	"$PosStrandFpkmNormalizedSums"
awk -F '\t' -v OFS='\t'	'{print $1, $2, $3, $4, $5, $6, $7 / ($3 - $2)}' \
		"$NegStrandMillionsNormalizedSums" >	"$NegStrandFpkmNormalizedSums"

## With the heavy lifting of normalization done, we can now easily
## take the maximum FPKM value for each isoform.
logr	"Re-Sorting Files for Merging"
sort -rnk7 "$PosStrandFpkmNormalizedSums" | sort -u -k4 | sort -k 1,1 -k2,2n > "$PosStrandFinalSorted"
sort -rnk7 "$NegStrandFpkmNormalizedSums" | sort -u -k4 | sort -k 1,1 -k2,2n > "$NegStrandFinalSorted"
wait

## Merge the two strands together. No sorting	is necessary yet, so save the bandwidth.
cat	"$PosStrandFinalSorted"	"$NegStrandFinalSorted" > "$StrandsMerged"

## Next, we use a supplementary python script to add common gene ID's
## so that we can filter out isoforms from our final file.
python "$ScriptDir"/convert_isoform.py -l "$ScriptDir"/refseq_to_common_id.txt \
			 -f "$StrandsMerged" \
			 -o "$FPKMCommonID"

## At long last, we sort in descending order and filter
sort -rnk7 "$FPKMCommonID" | sort -usk8 |
		sort -k1,1 -k2,2 > "$FinalFPKM"

awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6}' "$FinalFPKM" > "$FinalOut"

logr "Final Genes: ""$(wc -l "$FinalOut" | awk '{print $1}')"

## If we made it this far, we are done
logr "Done executing"
