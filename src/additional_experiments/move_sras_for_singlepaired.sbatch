#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH --mail-user=zama8258@colorado.edu

set -euxo pipefail

root_path=/scratch/Users/zama8258/virtual_spike_in/other_out/sras/
single_path=/scratch/Users/zama8258/virtual_spike_in/other_out/sras/single
paired_path=/scratch/Users/zama8258/virtual_spike_in/other_out/sras/paired
metadata=/scratch/Users/zama8258/virtual_spike_in/samples_with_pair.txt
single_samples=( $(<"$metadata"  grep single | awk '{print $1}') )
paired_samples=( $(<"$metadata"  grep paired | awk '{print $1}') )
for sample in "${single_samples[@]}"; do
		# rsync -P "$root_path"/"$sample".sra "$single_path"
		ln -s "$root_path"/"$sample".sra "$single_path"/"$sample".sra
done
for sample in "${paired_samples[@]}"; do
		# rsync -P "$root_path"/"$sample".sra "$paired_path"
		ln -s "$root_path"/"$sample".sra "$paired_path"/"$sample".sra
done

function finish {
		if command -v apprise &> /dev/null
		then
				apprise -t "Fiji Run Complete" -b "Move Data"
				exit
		fi
}
trap finish EXIT
