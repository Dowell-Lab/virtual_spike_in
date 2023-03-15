### compare_model_and_bulk.R ---
##
## Filename: compare_model_and_bulk.R
## Author: Zach Maas
## Created: Mon Sep 19 12:04:10 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## Compare the results from using bulk spikein normalization and the
## model to see where they differ.
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
library('DESeq2')
library('ggthemes')

## Load data
factors_bulk <- c(0.984, 0.81)
factors_model <- c(0.984, 0.81)
counts_file <- '~/Dropbox/phd/research/dna_lab/virtual_spike_in/pipeline_counts/counts/counts_merged.txt'
counts_ini <- read_delim(counts_file) %>%
    as_tibble()
counts <- counts_ini %>%
    transmute(control_1 = SRR5364303.sorted.sorted.bam,
              control_2 = SRR5364304.sorted.sorted.bam,
              treatment_1 = SRR5364305.sorted.sorted.bam,
              treatment_2 = SRR5364306.sorted.sorted.bam)
counts_matrix <- as.matrix(counts)
conditions <- factor(c("1", "1", "2", "2"))
coldata <- data.frame(row.names = colnames(counts_matrix), conditions)

## Run iterations
num_significant <- c()
res_tbls <- tibble()
run <- 0
for (factors in c(factors_bulk, factors_model)) {
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                  colData = coldata,
                                  design = ~conditions)
    res <- estimateSizeFactors(dds)
    sizeFactors(res) <- factors
    res <- results(DESeq(res))
    num_significant <- c(num_significant,
                         sum(res$padj < 0.001, na.rm=TRUE))
    res_tbl <- as_tibble(res) %>% mutate(factor = factor)
    res_tbls <- rbind(res_tbls, res_tbl)
}

significance_vs_factor <- tibble(factor = log2(factors),
                                 percent_significant = num_significant / 18093)
ggplot() +
    geom_line(data = significance_vs_factor,
              aes(x=factor,
                  y=percent_significant)) +
    geom_vline(aes(xintercept = 0), color="red") +
    ylim(0,1) +
    labs(title = "Proportion Significant vs Size Factor",
         x = "B / A Size Factor Ratio",
         y = "Proportion Significant (p < 0.01)") +
    theme_tufte()
plot_dir <- '/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/virtual_spike_in/out/new_plots/'
ggsave(paste0(plot_dir, "size_factor_plot.pdf"), width=8,height=8)

######################################################################
### compare_model_and_bulk.R ends here
