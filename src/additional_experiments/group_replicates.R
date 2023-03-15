### group_replicates.R --- Group replicates of interest
##
## Filename: group_replicates.R
## Author: Zach Maas
## Created: Mon Jan  2 15:30:08 2023 (-0700)
##
######################################################################
##
### Commentary:
##
## The replicates that we have are defined in the database by multiple
## identifiers, meaning that to group them we need to first aggregate
## identifiers and then group things from there.
##
######################################################################
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
##
######################################################################
##
### Code:

library('tidyverse')

baseDir <- '/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/'
dat <- read_delim(paste0('dat/rep_metadata.txt'))

## This accomplishes the first goal but not the second...
## Some experiments should have multiple sets of controls...
condition_agg <- dat %>%
    group_by(sample_name, paper_id) %>%
    summarise(condition=paste(condition_id, collapse="_")) %>%
    group_by(paper_id, condition) %>%
    summarise(samples=paste(sample_name, collapse="_")) %>%
    separate_rows(samples) %>%
    ungroup() %>%
    transmute(id=paste0(paper_id, "_", condition), samples=samples)
write_delim(condition_agg, paste0(baseDir, "dat/samples_grouped.txt"))

######################################################################
### group_replicates.R ends here
