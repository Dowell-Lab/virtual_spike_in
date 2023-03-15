### size_factor_titrations.R --- Size factor titration experiment
##
## Filename: size_factor_titrations.R
## Author: Zach Maas
## Created: Wed Sep 14 13:09:43 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## This file contains code to run a titration experiment using DESeq2
## size factors to show how differential expression results are highly
## sensitve to the size factors / normalization factors used. We do
## this by fixing the size factors for two paired replicates and
## varying the log-size factors for the other replicates from -4 to 4
## over a grid of size X, plotting a trajectory of how points change
## and become differential across the size factor space.
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
library('scico')

## Make parallel
library("BiocParallel")
register(MulticoreParam(8))

## The range to iterate over
log2fold_change_range <- 2
num_points <- 1##256
factors <- 2 ^ seq(-1.47,
                   1.04,
                   length.out = num_points)
## factors <- 2 ^ seq(-log2fold_change_range,
##                    log2fold_change_range,
##                    length.out = num_points)

## Load data
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

## Define DESeq2 Code
run_deseq_with_plot <- function(size_factors, plot_title) {
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                  colData = coldata,
                                  design = ~conditions)
    ## Use not isfalse here so that false is the key to use default
    ## size factors
    res <- estimateSizeFactors(dds)
    if (length(size_factors) == 4) {
        sizeFactors(res) <- size_factors
    }
    res <- results(DESeq(res))
    res_tbl <- as_tibble(res) %>% mutate(factor = factor,
                                         Geneid = counts_ini$Geneid,
                                         significant = padj < 0.01)

    ggplot(data = res_tbl) +
        geom_hline(yintercept=0, color = "grey", linetype="dashed") +
        geom_point(aes(x = baseMean, y = log2FoldChange), show.legend = FALSE) +
        scale_color_manual(values = c("black", "red")) +
        scale_x_log10(limits=c(1,50000)) +
        ylim(-8, 8) +
        labs(title=paste0(plot_title,
                          "\nNumber of Significant Points (p<0.01) = "),
             x="Count", y="Log2 Fold Change")
    return(res_tbl)
}

## Run iterations
plot_dir <- '/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/virtual_spike_in/out/new_plots/'
num_significant <- c()
res_tbls <- tibble()
run <- 0
for (factor in factors) {
    run <- run+1
    print(paste0("(", (run / num_points) * 100, "%) ",
                 "Model ", run, " of ", num_points))
    res <- run_deseq_with_plot(c(1, 1, factor, factor),
                               paste0("Size factor B / A = ", factor))
    res_tbls <- rbind(res_tbls, res)
    signif_count <- sum(res$padj < 0.01, na.rm=TRUE)
    num_significant <- c(num_significant,
                         signif_count)
    ggsave(paste0(plot_dir, "maplot_pdf/maplot_", run, ".pdf"), width=16,height=9)
    ggsave(paste0(plot_dir, "maplot_png/maplot_", run, ".png"), width=16,height=9)
}

## Make plots using collected data
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
         y = "Proportion Significant (p < 0.01)")
ggsave(paste0(plot_dir, "size_factor_plot.pdf"), width=8,height=8)

res_signif <- res_tbls %>% filter(significant == TRUE)
ggplot(data = res_signif) +
    geom_hline(yintercept=0, color = "grey", linetype="dashed") +
    geom_point(aes(x = baseMean, y = log2FoldChange,
                   color = log2(factor)), alpha=0.01) +
    scale_x_log10(limits=c(1,50000)) +
    ylim(-8, 8) +
    labs(title="Significance Front",
         x="Count", y="Log2 Fold Change") +
    guides(color=guide_legend(title="Log Normalization Factor"))
ggsave(paste0(plot_dir, "significance_front.jpg"), width=16,height=9)

## Run using the true values
## res_vsi <- run_deseq_with_plot(c(1, 0.890, 0.677, 1.095), "VSI Size Factors")
res_vsi <- run_deseq_with_plot(c(1, 0.890, 0.677, 1.095), "VSI Size Factors")
ggsave(paste0(plot_dir, "maplot_vsi.pdf"), width=16,height=9)

## Run using DESeq2 estimate
res_deseq <- run_deseq_with_plot(FALSE, "DESeq2 Size Factors")
ggsave(paste0(plot_dir, "maplot_deseq.pdf"), width=16,height=9)

## Show that significant values are a subset
res_deseq_signif <- res_deseq %>% subset(significant == TRUE)
res_vsi_signif <- res_vsi %>% subset(significant == TRUE)
res_deseq_overlap <- res_deseq %>%
    mutate(overlap = Geneid %in% res_vsi_signif$Geneid)
vsi_overlaps_deseq <- res_deseq %>%
    subset(Geneid %in% res_vsi_signif$Geneid)
proportion_vsi_preserved <- sum(vsi_overlaps_deseq$significant) / length(vsi_overlaps_deseq$significant)
## Subset plot... We want to show overlap dots and exceptions.

res_merge <- merge(res_vsi, res_deseq, by="Geneid") %>% as_tibble() %>%
    mutate(signif_merge = as_factor(significant.x + significant.y))
## Null is that p-values are not monotonic, so use spearman test with p<0.05
cor(res_merge$padj.x, res_merge$padj.y, method="spearman")
## Generate a diagnostic plot
plot_min <- -log10(min(res_merge$padj.x, res_merge$padj.y, na.rm=TRUE))
plot_max <- -log10(1)
ggplot(data = res_merge) +
    geom_point(aes(x=-log10(padj.y), y=-log10(padj.x), color=signif_merge)) +
    geom_vline(aes(xintercept=-log10(0.01)), color="red", linetype="dotted") +
    geom_hline(aes(yintercept=-log10(0.01)), color="red", linetype="dotted") +
    ## xlim(c(plot_max, plot_min)) +
    ## ylim(c(plot_max, plot_min)) +
    guides(color="none") +
    theme_minimal() +
    labs(x="DESeq2 -log10 p-adjusted", y="VSI -log10 p-adjusted",
         title="Genes called as significant with VSI and DESeq2 size factors (p < 0.01)")
ggsave(paste0(plot_dir, "significant_gene_comparison.pdf"), width=8,height=8)

ggplot(data = res_merge) +
    geom_hline(yintercept=0, color = "grey", linetype="dashed") +
    geom_point(aes(x = baseMean.y, y = log2FoldChange.y, color=signif_merge),
               show.legend = TRUE) +
    scale_color_manual(,
                       values = c("black", "red", "purple"),
                       labels = c("Not Significant", "DESeq2 Significant",
                                  "DESeq2 + VSI Significant")) +
    scale_x_log10(limits=c(1,50000)) +
    ylim(-8, 8) +
    theme_minimal() +
    labs(title=paste0("Differentially Expressed Genes in DESeq2 and VSI"),
         x="Count", y="Log2 Fold Change",
         color = "Significant Genes (p < 0.01)") +
    theme(legend.position="bottom")
ggsave(paste0(plot_dir, "deseq_comparison_plot.pdf"), width=16,height=9)

simulation_counts <- read_delim("/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/dat/deseq_variations/simulation_counts.txt",
                                col_names = c("Geneid", "count"))
res_count <- res_merge %>%
    left_join(simulation_counts, by="Geneid") %>%
    mutate(proportion=count/1000)

ggplot(data = res_count) +
    geom_hline(yintercept=0, color = "grey", linetype="dashed") +
    geom_point(aes(x = baseMean.y, y = log2FoldChange.y, color=proportion),
               show.legend = TRUE) +
    scale_color_gradient(low="#9E9E9E", high="#f15156", limits=c(0,1)) +
    scale_x_log10(limits=c(1,50000)) +
    theme_minimal() +
    ylim(-8, 8) +
    labs(title=paste0("Differentially Expressed Genes in DESeq2 Simulations"),
         x="Count", y="Log2 Fold Change",
         color = "Frequency where p<0.01") +
    theme(legend.position="bottom")
ggsave(paste0(plot_dir, "deseq_count_plot.pdf"), width=16,height=9)
ggsave(paste0(plot_dir, "deseq_count_plot.png"), width=16,height=9)

######################################################################
### size_factor_titrations.R ends here
