### random_deseq2_sample.R --- Run DESeq2 with sampled size factors
##
## Filename: random_deseq2_sample.R
## Author: Zach Maas
## Created: Wed Mar  1 14:29:08 2023 (-0700)
##
######################################################################
##
### Commentary:
##
## This script is intended to run many times in parallel. It provides
## a DESeq2 execution script that will run a single differential
## expression analysis using size factors sampled from a distribution
## estimated by the VSI model.
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
library('progress')

## Read in dataset as in compare_replicates.R
dm6_size <- 30000000 ## Transcriptome Depth
baseDir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/dat/"
res_vsi_human <- read_delim(paste0(baseDir, "all_virtual_human.txt"),
                            col_names = c("comparison", "factor_virtual_human", "var_virtual_human")
                            ) %>% distinct(comparison, .keep_all = TRUE)
res_linear_human <- read_delim(paste0(baseDir, "all_linear_human.txt"),
                               col_names = c("comparison", "factor_linear_human", "var_linear_human")
                               ) %>% distinct(comparison, .keep_all = TRUE)
res_vsi_drosophila <- read_delim(paste0(baseDir, "all_virtual_drosophila.txt"),
                                 col_names = c("comparison", "factor_virtual_drosophila", "var_virtual_drosophila")
                                 ) %>% distinct(comparison, .keep_all = TRUE)
res_linear_drosophila <- read_delim(paste0(baseDir, "all_linear_drosophila.txt"),
                                    col_names = c("comparison", "factor_linear_drosophila", "var_linear_drosophila")
                                    ) %>% distinct(comparison, .keep_all = TRUE)
millions_mapped <- read_delim(paste0(baseDir, "all_millionsmapped.txt"),
                              col_names = c("sample", "millionsmapped")
                              )
percent_mapped <- read_delim(paste0(baseDir, "../src/dukler/mapped_ratio.txt"),
                             col_names = c("sample", "percentage")
                             )
srr_mapping <- read_delim(paste0(baseDir, "experiment_to_srr.txt"),
                          col_names = c("experiment", "sample")
                          ) %>%
    separate("experiment", c("experiment"), extra="drop", sep="_") %>%
    distinct()
res_all <- res_vsi_human %>%
    full_join(res_linear_human, by = c("comparison")) %>%
    full_join(res_vsi_drosophila, by = c("comparison")) %>%
    full_join(res_linear_drosophila, by = c("comparison")) %>%
    mutate(comparison = str_remove_all(comparison, '.sorted')) %>%
    separate(comparison, c("sample", "reference"), "_vs_", FALSE) %>%
    inner_join(srr_mapping, by=c("sample"))

## Filter to our size factors
factor_1 <- log2((res_all %>% filter(sample == "SRR5364304"))$factor_virtual_human)
factor_2 <- log2((res_all %>% filter(sample == "SRR5364311"))$factor_virtual_human)
factor_3 <- log2((res_all %>% filter(sample == "SRR5364312"))$factor_virtual_human)
var_1 <- (res_all %>% filter(sample == "SRR5364304"))$var_virtual_human
var_2 <- (res_all %>% filter(sample == "SRR5364311"))$var_virtual_human
var_3 <- (res_all %>% filter(sample == "SRR5364312"))$var_virtual_human

## Read in DESeq2 count data
counts_file <- '~/Dropbox/phd/research/dna_lab/virtual_spike_in/pipeline_counts/counts/counts_merged.txt'
counts_ini <- read_delim(counts_file) %>%
    as_tibble()
counts <- counts_ini %>%
    transmute(control_1 = SRR5364303.sorted.sorted.bam,
              control_2 = SRR5364304.sorted.sorted.bam,
              treatment_1 = SRR5364311.sorted.sorted.bam,
              treatment_2 = SRR5364312.sorted.sorted.bam)
counts_matrix <- as.matrix(counts)
rownames(counts_matrix) <- counts_ini$Geneid
conditions <- factor(c("1", "1", "2", "2"))
coldata <- data.frame(row.names = colnames(counts_matrix), conditions)

## Set up for output
root_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/"
out_dir <- paste0(root_dir, "dat/deseq_variations/")

## For repeated runs
run_deseq <- function(size_factors) {
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData = counts_matrix,
                                                   colData = coldata,
                                                   design = ~conditions))
    ## Use not isfalse here so that false is the key to use default
    ## size factors
    ## res <- estimateSizeFactors(dds)
    ## if (length(size_factors) == 4) {
    sizeFactors(dds) <- size_factors
    ## }
    suppressMessages(res <- results(DESeq(dds)))
    res_tbl <- as_tibble(res) %>% mutate(Geneid = counts_ini$Geneid,
                                         significant = padj < 0.01)
    ## ggplot(data = res_tbl) +
    ##     geom_hline(yintercept=0, color = "grey", linetype="dashed") +
    ##     geom_point(aes(x = baseMean, y = log2FoldChange), show.legend = FALSE) +
    ##     scale_color_manual(values = c("black", "red")) +
    ##     scale_x_log10(limits=c(1,50000)) +
    ##     ylim(-8, 8) +
    ##     labs(title=paste0(plot_title,
    ##                       "\nNumber of Significant Points (p<0.01) = "),
    ##          x="Count", y="Log2 Fold Change")
    write_delim(res_tbl,
                paste0(out_dir,
                       "deseq_res_",
                       size_factors[2], "_",
                       size_factors[3], "_",
                       size_factors[4], ".txt"))
    return(res_tbl)
}
get_factors <- function() {
    factors <- c(1,
                 2^rnorm(1, factor_1, var_1),
                 2^rnorm(1, factor_2, var_2),
                 2^rnorm(1, factor_3, var_3))
    return(factors)
}

num_runs <- 1000
num_significant <- c()
pb <- progress_bar$new(total = num_runs,
                       format = "[:bar] :percent eta: :eta")
for (run in seq(num_runs)) {
    ## print(paste0("Run ", run))
    pb$tick()
    res <- run_deseq(get_factors())
    signif_count <- sum(res$padj < 0.01, na.rm=TRUE)
    num_significant <- c(num_significant,
     signif_count)
}

######################################################################
### random_deseq2_sample.R ends here
