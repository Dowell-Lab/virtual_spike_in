# fetch_data.sh --- run fastestq dump on data
#
# Filename: fetch_data.sh
# Author: Zach Maas
# Created: Thu Dec  8 09:43:28 2022 (-0700)
#
#

# Commentary:
#
#
#
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

script=/scratch/Users/zama8258/virtual_spike_in/src/dukler/fastestq_dump.sbatch
srrfile=/scratch/Users/zama8258/virtual_spike_in/other_out/sras/srrs.txt
outdir=/scratch/Users/zama8258/virtual_spike_in/other_out/sras
sbatch "$script" -f "$srrfile" -o "$outdir"

#
# fetch_data.sh ends here
