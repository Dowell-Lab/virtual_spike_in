#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/generic/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/generic/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --mail-user=zama8258@colorado.edu

# Standard bash settings
set -eo pipefail

while getopts ":ho:f:" opt; do
		case ${opt} in
				h )
						echo "Utility to Quickly and Reliably fetch SRA files."
						echo "Input: A text file with an SRR on each line"
						echo ""
						echo "Usage:"
						echo "  fastestq-dump -f srrs.txt -o /output/dir"
						echo "  - srrs.txt is a list of SRRs to fetch.  "
						echo "  - /output/dir is an existing directory  "
						echo "    to put your downloads into.           "
						exit 0
						;;
				f )
						echo "-f: ""$OPTARG"
						sraFile="$OPTARG"
						;;
				o )
						echo "-o: ""$OPTARG"
						outDir="$OPTARG"
						;;
				\? )
						echo "Invalid Option: -" "$OPTARG" 1>&2
						exit 1
						;;
		esac
done
shift $((OPTIND -1))

echo "Working in: ""$PWD"
if [ ! -f "$sraFile" ]; then
		echo "File -f does not exist."
		exit 1
fi
if [ ! -d "$outDir" ]; then
		echo "Directory -o does not exist."
		exit 1
fi

# Load sra toolkit for fastq-dump
module load sra

# Temporary files
ariaList=$(mktemp /tmp/fastestq-dump.XXXXXX)

# Parse the list of SRRs
parallel --trim r -a "$sraFile" echo -e '$('srapath {}')'"'\'n'\'tout={}.sra" > "$ariaList"
cat "$ariaList"

# The vital features here are max-concurrent-downloads and
# max-connection-per-server, which both limit the number of simultaneous
# connections that we attempt to make. Here I've limited it to 8 (equivalent to
# running fasterq-dump or parallel-fastq-dump with 8 threads). You can change
# this, just be sure to be considerate.
num_max_dl=8

# You'll need to point the script at your aria2 executable. Mine is in my path,
# so it's not an issue. Otherwise, you need to change this variable to point to
# the location of aria2 on FIJI
aria2bin=$(which aria2c)

# Download so fast.
"$aria2bin" \
		--dir="$outDir"  \
		--input-file="$ariaList" \
		--continue=true \
		--enable-http-pipelining=true \
		--max-concurrent-downloads="$num_max_dl" \
		--connect-timeout=60  \
		--max-connection-per-server="$num_max_dl" \
		--split=4 \
		--file-allocation=prealloc \
		--disk-cache=512M \
		--min-split-size=5M  \
		--human-readable=true \
		--enable-color=false
