### compare_experiments.R --- Compare experimental output
##
## Filename: compare_experiments.R
## Author: Zach Maas
## Created: Fri Dec  2 12:35:16 2022 (-0700)
##
######################################################################
##
### Commentary:
##
## This file contains code to analyze the results of several
## additional experiments not present in the original analysis for
## this project.
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

library("tidyverse")
library("ggthemes")
library('patchwork')

## Read in data
dm6_size <- 143726002
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
    inner_join(millions_mapped, by = c("sample")) %>%
    inner_join(millions_mapped, by = c("reference" = "sample")) %>%
    rename(c(depth_sample = millionsmapped.x,
             depth_reference = millionsmapped.y)) %>%
    inner_join(srr_mapping, by=c("sample")) %>%
    mutate(factor_mapped = depth_reference / depth_sample,
		proportion_reference = depth_reference / dm6_size,
		proportion_sample = depth_sample / dm6_size)
nrow(res_all)
n_expt <- length(unique(srr_mapping$experiment))
n_sample <- length(unique(res_all$sample))

## Generate comparison plots within methods
plt_1 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_human),
        y = log2(factor_virtual_human),
        color = experiment,
				alpha = depth_sample,
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "VSI/Linear (Human)",
        x = "Log2(Size Factor Linear)",
        y = "Log2(Size Factor VSI)",
        color = "Experiment"
    )
ggsave("linear_vs_virt_human.pdf", width = 9, height = 9)

plt_2 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_drosophila),
        y = log2(factor_virtual_drosophila),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "VSI/Linear (Drosophila)",
        x = "Log2(Size Factor Linear)",
        y = "Log2(Size Factor VSI)",
        color = "Experiment"
    )
ggsave("linear_vs_virt_drosophila.pdf", width = 9, height = 9)

plt_3 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_drosophila),
        y = log2(factor_mapped),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Mapped/Linear Exogenous",
        x = "Log2(Size Factor Linear)",
        y = "Log2(Size Factor Mapped)",
        color = "Experiment"
    )
ggsave("mapped_vs_linear_drosophila.pdf", width = 9, height = 9)

plt_4 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_virtual_drosophila),
        y = log2(factor_mapped),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Mapped/VSI Exogenous",
        x = "Log2(Size Factor VSI)",
        y = "Log2(Size Factor Mapped)",
        color = "Experiment"
    )
ggsave("mapped_vs_virt_drosophila.pdf", width = 9, height = 9)

plt_5 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_human),
        y = log2(factor_mapped),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Mapped/Linear 3'",
        x = "Log2(Size Factor Linear)",
        y = "Log2(Size Factor Mapped)",
        color = "Experiment"
    )
ggsave("mapped_vs_linear_human.pdf", width = 9, height = 9)

plt_6 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_mapped),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Mapped/VSI 3'",
        x = "Log2(Size Factor VSI)",
        y = "Log2(Size Factor Mapped)",
        color = "Experiment"
    )
ggsave("mapped_vs_virt_human.pdf", width = 9, height = 9)

## Generate plots comparing choice of regions
plt_7 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_human),
        y = log2(factor_linear_drosophila),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Exogenous/3' (Linear)",
        x = "Log2(Size Factor 3')",
        y = "Log2(Size Factor Drosophila)",
        color = "Experiment"
    )
ggsave("human_vs_drosophila_linear.pdf", width = 9, height = 9)

plt_8 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_virtual_drosophila),
        color = experiment,
				alpha = depth_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed() +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Exogenous/3' (VSI)",
        x = "Log2(Size Factor 3')",
        y = "Log2(Size Factor Drosophila)",
        color = "Experiment"
    )
ggsave("human_vs_drosophila_virt.pdf", width = 9, height = 9)

plt_1 + plt_2 + plt_3 + plt_4 + plt_5 + plt_6 + plt_7 + plt_8 +
    plot_layout(guides = "collect", ncol=4) & theme(legend.position = "right")
ggsave("comparisons_merged.pdf", width = 16, height = 9)

print(cor(res_all$factor_linear_human, res_all$factor_linear_drosophila, use = "pairwise"))
print(cor(res_all$factor_virtual_human, res_all$factor_virtual_drosophila, use = "pairwise"))

ggplot(data = res_all) +
		geom_point(aes(x=depth_reference,
                   y=depth_sample,
                   color=experiment)) +
    geom_vline(xintercept=1307030, color="red", linetype="dashed") +
    geom_hline(yintercept=1307030, color="red", linetype="dashed") +
    geom_abline(slope=1, intercept=0, color="black", linetype="dotted") +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    labs(
        title = "Published Spikeins Have High Variability",
        x = "Depth Reference",
        y = "Depth Sample",
        color="Experiment"
    )
ggsave("spikein_consistency.pdf", width = 9, height = 9)

ggplot(data = res_all) +
		geom_point(aes(x=pmin(depth_sample, depth_reference),
                   y=abs(log2(factor_virtual_human)-
                         log2(factor_virtual_drosophila)),
                   color=experiment)) +
    geom_vline(xintercept=1307030, color="red", linetype="dotted") +
    geom_vline(xintercept=3485415, color="blue", linetype="dotted") +
    theme_tufte() +
    theme(text=element_text(family="Helvetica")) +
    scale_x_log10() +
    labs(
        title = "Sample Consistency Varies with Spikein Depth",
        x = "Spikein Depth",
        y = "Absolute Difference Between Spikein and 3'",
        color="Experiment"
    )
ggsave("estimate_variability.pdf", width = 9, height = 9)

######################################################################
### compare_experiments.R ends here
