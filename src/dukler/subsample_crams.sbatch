#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=32gb
#SBATCH --mail-user=zama8258@colorado.edu

# subsample_crams.sbatch --- Subsample crams for tiration experiment
#
# Filename: subsample_crams.sbatch
# Author: Zach Maas
# Created: Mon Apr 18 12:29:12 2022 (-0600)
#
#

# Commentary:
#
#
# This script subsamples a set of crams of interest (3 WT samples from
# an experiment with high depth) in order to establish a best case on
# the kinds of spikein uncertainty we can expect as we decrease the
# depth of the spikein. We do n-many subsamples to get a statistical
# bound on our subsampling. Arbitrarily chosen as n=128
# 
# 
# this program is free software: you can redistribute it and/or modify
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

set -euxo pipefail

module load samtools

FASTA=/scratch/Shares/dowell/genomes/dm6/dm6.fa
BASEDIR=/scratch/Users/zama8258/virtual_spike_in/other_out/subsample_experiment
NP=32
# parallel -j "$NP" --progress samtools view -C -T "$FASTA" -s {3}{1} {2} '>' {2.}_{1}_{3}.cram ::: \
		#  				 .5 .25 .125 .0625 .03125 .015625 ::: \
		#  				 "$BASEDIR"/crams/*.cram ::: \
		# 				 $(seq 128)
# parallel -j "$NP" --progress samtools index {} ::: \
		#  				 "$BASEDIR"/crams/*.cram
parallel -j "$NP" --progress samtools flagstat {} '>' {}.stats ::: \
				 "$BASEDIR"/crams/*.cram
tail -n +0 *.stats > all_mapstats.stat

#
# subsample_crams.sbatch ends here
