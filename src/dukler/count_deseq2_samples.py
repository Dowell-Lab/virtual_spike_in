# count_deseq2_samples.py --- Count variations from our DESeq2 samples
#
# Filename: count_deseq2_samples.py
# Author: Zach Maas
# Created: Thu Mar  2 14:11:10 2023 (-0700)
#
#

# Commentary:
#
#
# This script uses Python to count the number of times a gene is
# called as significant by DESeq2 in our random sampling process
# across the space of feasible normalization parameters. This is done
# in Python because the data structures are easier to work with.
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

from collections import Counter
from glob import glob
from tqdm import tqdm

root_dir = "/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in"
data_dir = f"{root_dir}/dat/deseq_variations"

res_file_list = glob(f"{data_dir}/deseq_res_*.txt")
count = Counter()

for res in tqdm(res_file_list):
    with open(res, "r") as curr_file:
        curr_file.readline()  # Skip header
        for line in curr_file:
            line_list = line.strip().split()
            try:
                if float(line_list[5]) < 0.01:
                    count.update([line_list[6]])
            except:
                continue

with open(f"{data_dir}/simulation_counts.txt", "w") as out_file:
    for k, count in count.most_common():
        out_file.write(f"{k}\t{count}\n")

#
# count_deseq2_samples.py ends here
