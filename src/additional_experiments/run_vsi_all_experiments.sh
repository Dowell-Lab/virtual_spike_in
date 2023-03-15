#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1gb
#SBATCH --mail-user=zama8258@colorado.edu
# run_vsi_all_experiments.sh --- All in one script
#
# Filename: run_vsi_all_experiments.sh
# Author: Zach Maas
# Created: Fri Nov 18 10:12:23 2022 (-0700)
#
#

# Commentary:
#
#
# This script serves to move data from multiple previously published
# experiments totaling 288 samples, and then process it using the
# virtual spike-in script. We do this for both the drosophila data
# (separately analyzed) and the human data (from the database).
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
#
#

# Code:

set -uxo pipefail

# Experiments to analyze
experiments=(
		Aoi2020nelf
		Barbieri2020rapid
		Birkenheuer2018herpes
		Birkenheuer2020rna
		Dukler2017nascent
		Fan2020drb
		Jaeger2020selective
		Leroy2019ledgf
		Liu2021transcription
		Rao2017cohesin
		Santoriello2020rna
		Sendinc2019pcif1
		Takahashi2020role
		Vihervaara2021stress
)

# Set up directory structure
scratch="/scratch/Users/zama8258"
baseDir="$scratch"/vsi_additional_experiments
mkdir -p "$baseDir"
# Make individual experiment folders
for experiment in "${experiments[@]}"; do
		mkdir -p "$baseDir"/"$experiment"_human
		mkdir -p "$baseDir"/"$experiment"_drosophila
done

# Copy crams to directories
# TODO Check this directory
drosophila_cram_dir="$scratch"/virtual_spike_in/other_out/nascentflow_out/mapped/crams
for experiment in "${experiments[@]}"; do
		echo "Copying Human"
		rsync -Pra \
					"/Shares/dbnascent/""$experiment""/crams/" \
					"$baseDir"/"$experiment"_human
		expt_srrs=(
				$(find "$baseDir"/"$experiment"_human -type f | grep .cram$ | xargs -I % basename % .sorted.cram)
		)
		echo "Copying Drosophila"
		for srr in "${expt_srrs[@]}"; do
				rsync -Pa \
							"$drosophila_cram_dir"/"$srr"".sorted.cram" \
							"$drosophila_cram_dir"/"$srr"".sorted.cram.crai" \
							"$baseDir"/"$experiment"_drosophila
		done
done

#
# run_vsi_all_experiments.sh ends here
