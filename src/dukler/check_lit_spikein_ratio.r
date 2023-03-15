### check_lit_spikein_ratio.r --- Check the ratio of mapped reads
##
## Filename: check_lit_spikein_ratio.r
## Author: Zach Maas
## Created: Wed Apr 13 09:26:37 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## We're interested in knowing how the proportion of reads spiked-in
## varies across different experiments, so we took 11 different
## experiments that do a drosophila spike-in, mapped to dm6, and
## looked at the percentage of reads that hisat2 mapped.
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
library('ggplot2')
library('ggthemes')
library('extrafont')

baseDir <- '/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in'
df_percent <- read_delim(paste0(baseDir, '/src/dukler/mapped_ratio.txt'),
                         col_names = c('run', 'percentage'))

dm6_size <- 30000000
df_mapped <- read_delim(paste0(baseDir, '/dat/all_millionsmapped.txt'),
                        col_names = c('run', 'millionsmapped')) %>%
    mutate(coverage = (millionsmapped * 150) / dm6_size)

srr_mapping <- read_delim(paste0(baseDir, "/dat/replicates/experiment_to_srr.txt"),
                          col_names = c("experiment", "sample")
                          ) %>%
    separate("experiment", c("experiment"), extra="drop", sep="_") %>%
    distinct()

df_joined <- left_join(df_mapped, df_percent, by="run") %>%
    na.omit()  %>%
    filter(run %in% srr_mapping$sample) %>%
    mutate(cov_comparison = coverage / percentage)

n_expt <- length(unique(srr_mapping$experiment))
n_sample <- length(unique(df_joined$run))

plotDir <- paste0(baseDir, '/dukler_out/plots/')
ggplot(data=df_percent) +
    geom_histogram(aes(x=percentage)) +
    labs(x = 'Percent Mapped to Spike-In',
         y = 'Count',
         title = paste0('Spike-In Mapping (',
                        n_expt, " Publications, n=", n_sample, ")"))
print(plotDir)
ggsave(paste0(plotDir, "spikein_proportions.pdf"))

ggplot(data=df_joined) +
    geom_histogram(aes(x=millionsmapped)) +
    labs(x = 'Reads Mapped to Spike-In (Log Scale)',
         y = 'Count',
         title = paste0('Spike-In Mapping (',
                        n_expt, " Publications, n=", n_sample, ")")) +
    scale_x_log10(breaks=c(1e4,1e5,1e6,1e7,1e8))
print(plotDir)
ggsave(paste0(plotDir, "spikein_millionsmapped.pdf"))

ggplot(data=df_joined) +
    geom_density(aes(x=coverage)) +
    ## Assume 30mil deep experiment on human hg38 genome:
    geom_vline(aes(xintercept=10, color="red")) +
    labs(x = 'Average Coverage per Genome Base',
         y = 'Count',
         title = paste0('Spike-In Genome Coverage (',
                        n_expt, " Publications, n=", n_sample, ")")) +
    theme(text=element_text(family="Helvetica"), legend.position="none")
print(plotDir)
ggsave(paste0(plotDir, "spikein_millionsmapped_coverage.pdf"), width=8, height=6)

ggplot(data=df_joined) +
    geom_point(aes(x=percentage, y=coverage)) +
    labs(x = 'Percent of Reads Mapped to Spike-In',
         y = 'Average Coverage per Genome Base',
         title = paste0('Average Coverage vs Mapping Percentage (',
                        n_expt, " Publications, n=", n_sample, ")")) +
    theme(legend.position="bottom", legend.justification="right",
          text=element_text(family="Helvetica"))
print(plotDir)
ggsave(paste0(plotDir, "spikein_coverage_vs_percentage.pdf"), width=8, height=6)

######################################################################
### check_lit_spikein_ratio.r ends here
