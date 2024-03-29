# run.sh --- Run Nextflow Pipeline
#
# Filename: run.sh
# Author: Zach Maas
# Created: Wed Apr  6 12:11:09 2022 (-0600)
#
#

# Commentary:
#
#
# Run the nextflow pipeline on test data
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

eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv activate virtual_spike_in
nextflow run main.nf -profile example -resume

#
# run.sh ends here
