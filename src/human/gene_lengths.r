### gene_lengths.r --- Plot showing gene lengths
##
## Filename: gene_lengths.r
## Author:
## Created: Mon Jul 26 15:08:57 2021 (-0600)
##
######################################################################
##
### Commentary:
##
## This file contains code to show how the number of genes available
## for the virtual-spike-in method drops with increasing timepoint:w
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

## Import and parse
path <- '/home/zach/Dropbox/phd/research/dna_lab/virtual_spike_in/dat/hg38_refseq.bed'
refseq <- read_delim(path, delim='\t', col_names=FALSE) %>%
    mutate(length = X3 - X2, hours = length / (2000*60))

## TODO This needs to use isoform filtered refseq, maybe?
for (hour in c(1, 2, 3, 4, 5, 6)) {
    print(nrow(filter(refseq, hours >= hour)))
}

ggplot(refseq, aes(x=hours)) + stat_ecdf(geom = "step") +
    scale_y_reverse() +
    scale_y_log10() +
    xlim(0,5) +
    labs(title="Proportion of Remaining Usable Spike-In Genes Over Time",
         x="Hours Since Treatment", y="Remaining Proportion of Genes")
ggsave('/home/zach/Dropbox/phd/research/dna_lab/virtual_spike_in/figs/length_distr.pdf')

refseq %>% filter(length > 3000*160) %>% distinct(X2, .keep_all=TRUE) %>% distinct(X3, .keep_all=TRUE)

######################################################################
### gene_lengths.r ends here
