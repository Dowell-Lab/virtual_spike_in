### compare_replicates.R --- Compare experimental output
##
## Filename: compare_replicates.R
## Author: Zach Maas
## Created: Fri Dec  2 12:35:16 2022 (-0700)
##
######################################################################
##
### Commentary:
##
## This file contains code to analyze the results of several
## additional replicates not present in the original analysis for
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
library('lubridate')

## Read in data
## dm6_size <- 143726002 ## Whole Genome Depth
dm6_size <- 30000000 ## Transcriptome Depth
experiments <- c("Aoi2020nelf", "Barbieri2020rapid",
                 "Birkenheuer2018herpes", "Birkenheuer2020rna",
                 "Dukler2017nascent", "Fan2020drb",
                 "Jaeger2020selective", "Leroy2019ledgf",
                 "Liu2021transcription", "Rao2017cohesin",
                 "Santoriello2020rna", "Sendinc2019pcif1",
                 "Takahashi2020role", "Vihervaara2021stress")
## baseDir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in/dat/replicates/"
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
time_mapping <- read_delim(paste0(baseDir, "timepoint_metadata.txt"),
                           col_names = c("sample", "duration", "unit")) %>%
    transmute(sample = sample,
              time = coalesce(as.duration(paste0(duration, unit)),
                              duration("0m"))) %>%
    distinct()
res_all <- res_vsi_human %>%
    full_join(res_linear_human, by = c("comparison")) %>%
    full_join(res_vsi_drosophila, by = c("comparison")) %>%
    full_join(res_linear_drosophila, by = c("comparison")) %>%
    mutate(comparison = str_remove_all(comparison, '.sorted')) %>%
    separate(comparison, c("sample", "reference"), "_vs_", FALSE) %>%
    left_join(millions_mapped, by = c("sample")) %>%
    left_join(millions_mapped, by = c("reference" = "sample")) %>%
    left_join(percent_mapped, by = c("sample")) %>%
    left_join(percent_mapped, by = c("reference" = "sample")) %>%
    left_join(time_mapping, by = c("sample")) %>%
    left_join(time_mapping, by = c("reference" = "sample")) %>%
    rename(c(depth_sample = millionsmapped.x,
             depth_reference = millionsmapped.y,
             percent_sample = percentage.x,
             percent_reference = percentage.y,
             time_sample = time.x,
             time_reference = time.y)) %>%
    inner_join(srr_mapping, by=c("sample")) %>%
    mutate(factor_mapped = depth_reference / depth_sample,
           proportion_reference = depth_reference / dm6_size,
           proportion_sample = depth_sample / dm6_size,
           coverage_sample = proportion_sample*150,
           coverage_reference = proportion_reference*150) %>%
    filter(experiment %in% experiments)
conds_to_factor <- function(good_timepoint, good_depth) {
    if (good_timepoint == 1 & good_depth == 1) {
        return("All Assumptions Met")
    } else if (good_timepoint == 1 & good_depth == 0) {
        return("Low Depth")
    } else if (good_timepoint == 0 & good_depth == 1) {
        return("Long Timepoint")
    } else {
        return ("All Assumptions Violated")
    }
}
res_all <- res_all %>% mutate(good_timepoint =
                                  coalesce(as.integer(
                                      time_sample < duration("1h")), 0),
                              good_depth = coalesce(as.integer(
                                  coverage_sample > 10), 0),
                              both = as.factor(
                                  mapply(conds_to_factor,
                                         as.integer(good_timepoint),
                                         as.integer(good_depth))))
nrow(res_all)
n_expt <- length(unique(res_all$experiment))
n_sample <- length(unique(res_all$sample))

## Generate comparison plots within methods
plt_1 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_human),
        y = log2(factor_virtual_human),
        color = experiment,
        alpha = coverage_sample,
        )) +
    geom_errorbar(aes(x=log2(factor_linear_human),
                      ymin=log2(factor_virtual_human)-var_virtual_human,
                      ymax=log2(factor_virtual_human)+var_virtual_human,
                      color=experiment,
                      alpha=coverage_sample)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Method Similarity (3'): VSI/Linear",
        x = "Log2(Size Factor Linear)",
        y = "Log2(Size Factor VSI)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("linear_vs_virt_human.pdf", width = 9, height = 9)

plt_2 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_drosophila),
        y = log2(factor_virtual_drosophila),
        color = experiment,
        alpha = coverage_sample
    )) +
    geom_errorbar(aes(x=log2(factor_linear_drosophila),
                      ymin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                      ymax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                      color=experiment,
                      alpha=coverage_sample)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Method Similarity (Drosophila): VSI/Linear",
        x = "Log2(Size Factor Linear)",
        y = "Log2(Size Factor VSI)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("linear_vs_virt_drosophila.pdf", width = 9, height = 9)

plt_3 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_drosophila),
        y = log2(factor_mapped),
        color = experiment,
        alpha = coverage_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Exogenous Estimate Similarity (Linear/Point)",
        x = "Log2(Size Factor Linear Model)",
        y = "Log2(Size Factor Point Estimate)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("mapped_vs_linear_drosophila.pdf", width = 9, height = 9)

plt_4 <- ggplot(data = res_all %>% filter(both == "All Assumptions Met")) +
    geom_point(aes(
        x = log2(factor_virtual_drosophila),
        y = log2(factor_mapped),
        color = experiment,
        alpha = coverage_sample
    )) +
    geom_errorbarh(aes(y=log2(factor_mapped),
                       xmin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                       xmax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                       color=experiment,
                       alpha=coverage_sample)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_cartesian(xlim=c(-5,5), ylim=c(-2.5,2.5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Exogenous Estimate Similarity (VSI/Point)",
        x = "Log2(Size Factor VSI Exogenous)",
        y = "Log2(Size Factor Exogenous Point Estimate)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("mapped_vs_virt_drosophila.pdf", width = 16, height = 9)

plt_5 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_human),
        y = log2(factor_mapped),
        color = experiment,
        alpha = coverage_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Mapped Ratio / 3' Linear",
        x = "Log2(Size Factor Linear Model)",
        y = "Log2(Size Factor Point Estimate)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("mapped_vs_linear_human.pdf", width = 9, height = 9)

plt_6 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_mapped),
        color = experiment,
        alpha = coverage_sample
    )) +
    geom_errorbarh(aes(y=log2(factor_mapped),
                       xmin=log2(factor_virtual_human)-var_virtual_human,
                       xmax=log2(factor_virtual_human)+var_virtual_human,
                       color=experiment,
                       alpha=coverage_sample)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Mapped Ratio / 3' VSI",
        x = "Log2(Size Factor VSI 3')",
        y = "Log2(Size Factor Point Estimate)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("mapped_vs_virt_human.pdf", width = 9, height = 9)

## Generate plots comparing choice of regions
plt_7 <- ggplot(data = res_all) +
    geom_point(aes(
        x = log2(factor_linear_human),
        y = log2(factor_linear_drosophila),
        color = experiment,
        alpha = coverage_sample
    )) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Linear 3' / Linear Drosophila",
        x = "Log2(Size Factor Linear 3')",
        y = "Log2(Size Factor Linear Drosophila)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("human_vs_drosophila_linear.pdf", width = 9, height = 9)

plt_8 <- ggplot(data = res_all %>% filter(both == "All Assumptions Met")) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_virtual_drosophila),
        color = experiment,
        )) +
    geom_errorbar(aes(x=log2(factor_virtual_human),
                      ymin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                      ymax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                      color=experiment)) +
    geom_errorbarh(aes(y=log2(factor_virtual_drosophila),
                       xmin=log2(factor_virtual_human)-var_virtual_human,
                       xmax=log2(factor_virtual_human)+var_virtual_human,
                       color=experiment)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "VSI 3' / VSI Drosophila",
        x = "Log2(Size Factor VSI 3')",
        y = "Log2(Size Factor VSI Drosophila)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("human_vs_drosophila_virt.pdf", width = 9, height = 9)

ggplot(data = res_all %>% filter(both != "All Assumptions Met")) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_virtual_drosophila),
        color = experiment,
        )) +
    geom_errorbar(aes(x=log2(factor_virtual_human),
                      ymin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                      ymax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                      color=experiment)) +
    geom_errorbarh(aes(y=log2(factor_virtual_drosophila),
                       xmin=log2(factor_virtual_human)-var_virtual_human,
                       xmax=log2(factor_virtual_human)+var_virtual_human,
                       color=experiment)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "VSI 3' / VSI Drosophila",
        x = "Log2(Size Factor VSI 3')",
        y = "Log2(Size Factor VSI Drosophila)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("human_vs_drosophila_virt_invert.pdf", width = 9, height = 9)
ggplot(data = res_all %>%
           filter(coverage_sample <= 10)) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_virtual_drosophila),
        color = experiment,
        )) +
    geom_errorbar(aes(x=log2(factor_virtual_human),
                      ymin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                      ymax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                      color=experiment)) +
    geom_errorbarh(aes(y=log2(factor_virtual_drosophila),
                       xmin=log2(factor_virtual_human)-var_virtual_human,
                       xmax=log2(factor_virtual_human)+var_virtual_human,
                       color=experiment)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "VSI 3' / VSI Drosophila",
        x = "Log2(Size Factor VSI 3')",
        y = "Log2(Size Factor VSI Drosophila)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("human_vs_drosophila_virt_invert_low.pdf", width = 9, height = 9)
ggplot(data = res_all %>%
           filter(time_sample > duration("1h"))) +
    geom_point(aes(
        x = log2(factor_virtual_human),
        y = log2(factor_virtual_drosophila),
        color = experiment,
        )) +
    geom_errorbar(aes(x=log2(factor_virtual_human),
                      ymin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                      ymax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                      color=experiment)) +
    geom_errorbarh(aes(y=log2(factor_virtual_drosophila),
                       xmin=log2(factor_virtual_human)-var_virtual_human,
                       xmax=log2(factor_virtual_human)+var_virtual_human,
                       color=experiment)) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
    geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
    geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
    coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "VSI 3' / VSI Drosophila",
        x = "Log2(Size Factor VSI 3')",
        y = "Log2(Size Factor VSI Drosophila)",
        color = "Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("human_vs_drosophila_virt_invert_long.pdf", width = 9, height = 9)

plt_1 + plt_2 + plt_3 + plt_4 + plt_5 + plt_6 + plt_7 + plt_8 +
    plot_layout(guides = "collect", ncol=2) & theme(legend.position = "right")
ggsave("comparisons_merged.pdf", width = 10, height = 16)

plt_4 + plt_8 + plot_layout(guides = "collect", ncol=2) & theme(legend.position = "right")
ggsave("secondary_comparison.pdf", width = 16, height = 9)

print(cor(res_all$factor_linear_human, res_all$factor_linear_drosophila, use = "pairwise"))
print(cor(res_all$factor_virtual_human, res_all$factor_virtual_drosophila, use = "pairwise"))

################################
## plot based on data quality ##
################################
plot_with_assumption_colors <- function(color_values, levels, plot_file, label, error_set) {
    if (label) {
        labels <- labs(
            title = "VSI 3' / VSI Drosophila - goodness",
            x = "Log2(Size Factor VSI 3')",
            y = "Log2(Size Factor VSI Drosophila)",
            color = "Sample Goodness",
            alpha = "Sample Goodness"
        )
        guides <- guides()
        size <- 2
    } else {
        labels <- labs(title="", x="", y="")
        guides <- guides(color="none", alpha="none")
        size <- 4
    }
    res_sorted <- res_all %>%
        mutate(both = factor(both, levels=levels)) %>%
        arrange(both)
    ## Universal Plotting
    plt <- ggplot(data = res_sorted) +
        geom_point(aes(
            x = log2(factor_virtual_human),
            y = log2(factor_virtual_drosophila),
            color = both), size=size)
    ## Conditional Error Bars
    if (length(error_set) > 0) {
        res_filt <- res_sorted %>% filter(both %in% error_set)
        plt <- plt +
            geom_errorbar(data = res_filt,
                          aes(x=log2(factor_virtual_human),
                              ymin=log2(factor_virtual_drosophila)-var_virtual_drosophila,
                              ymax=log2(factor_virtual_drosophila)+var_virtual_drosophila,
                              color=both)) +
            geom_errorbarh(data = res_filt,
                           aes(y=log2(factor_virtual_drosophila),
                               xmin=log2(factor_virtual_human)-var_virtual_human,
                               xmax=log2(factor_virtual_human)+var_virtual_human,
                               color=both))
    }
    ## Add decoration
    plt + scale_color_manual(values = color_values) +
        geom_vline(xintercept = 0, linetype = 'dotted', color = "black") +
        geom_hline(yintercept = 0, linetype = 'dotted', color = "black") +
        geom_abline(slope = 1, intercept=0, linetype = 'dashed', color = "red") +
        coord_fixed(xlim=c(-5,5), ylim=c(-5,5)) +
        theme_tufte() +
        labels + guides +
            theme(text=element_text(family="Helvetica", size=16),
                  legend.position="top")
    ggsave(plot_file, width = 9, height = 9)
    return(plt)
}
##############################
## Plots without error bars ##
##############################
g1 <- plot_with_assumption_colors(c("All Assumptions Met" = "darkgreen",
                                    "All Assumptions Violated" = "red",
                              "Low Depth" = "orange",
                              "Long Timepoint" = "blue"),
                            c("All Assumptions Met",
                              "All Assumptions Violated",
                              "Low Depth",
                              "Long Timepoint"),
                            "greyout_full.pdf",
                            TRUE,
                            "")
g2 <- plot_with_assumption_colors(c("All Assumptions Met" = "darkgreen",
                                    "All Assumptions Violated" = "aliceblue",
                                    "Low Depth" = "aliceblue",
                                    "Long Timepoint" = "aliceblue"),
                                  c("Long Timepoint",
                              "All Assumptions Violated",
                              "Low Depth",
                              "All Assumptions Met"),
                            "greyout_all.pdf",
                            FALSE,
                            "")
g3 <- plot_with_assumption_colors(c("All Assumptions Met" = "aliceblue",
                                    "All Assumptions Violated" = "red",
                                    "Low Depth" = "aliceblue",
                                    "Long Timepoint" = "aliceblue"),
                                  c("All Assumptions Met",
                              "Long Timepoint",
                              "Low Depth",
                              "All Assumptions Violated"),
                            "greyout_none.pdf",
                            FALSE,
                            "")
g4 <- plot_with_assumption_colors(c("All Assumptions Met" = "aliceblue",
                                    "All Assumptions Violated" = "aliceblue",
                                    "Low Depth" = "orange",
                                    "Long Timepoint" = "aliceblue"),
                                  c("All Assumptions Met",
                              "All Assumptions Violated",
                              "Long Timepoint",
                              "Low Depth"),
                            "greyout_low.pdf",
                            FALSE,
                            "")
g5 <- plot_with_assumption_colors(c("All Assumptions Met" = "aliceblue",
                                    "All Assumptions Violated" = "aliceblue",
                                    "Low Depth" = "aliceblue",
                                    "Long Timepoint" = "blue"),
                            c("All Assumptions Met",
                              "All Assumptions Violated",
                              "Low Depth",
                              "Long Timepoint"),
                            "greyout_long.pdf",
                            FALSE,
                            "")
## layout <- "
## ##AAAA
## ##AAAA
## BBAAAA
## BBAAAA
## CCDDEE
## CCDDEE
## "
## g1 + g2 + g3 + g4 + g5 +
##     plot_layout(design = layout, guides="collect")
## ggsave("merged_assumptions.pdf", width=16, height=9)

############################
## Replot with error bars ##
############################
e1 <- plot_with_assumption_colors(c("All Assumptions Met" = "darkgreen",
                                    "All Assumptions Violated" = "red",
                                    "Low Depth" = "orange",
                                    "Long Timepoint" = "blue"),
                                  c("All Assumptions Met",
                                    "All Assumptions Violated",
                                    "Low Depth",
                                    "Long Timepoint"),
                                  "err_greyout_full.pdf",
                                  TRUE,
                                  c("All Assumptions Met",
                                    "All Assumptions Violated",
                                    "Low Depth",
                                    "Long Timepoint"))
e2 <- plot_with_assumption_colors(c("All Assumptions Met" = "darkgreen",
                                    "All Assumptions Violated" = "aliceblue",
                                    "Low Depth" = "aliceblue",
                                    "Long Timepoint" = "aliceblue"),
                                  c("Long Timepoint",
                                    "All Assumptions Violated",
                                    "Low Depth",
                                    "All Assumptions Met"),
                                  "err_greyout_all.pdf",
                                  FALSE,
                                  c("All Assumptions Met"))
e3 <- plot_with_assumption_colors(c("All Assumptions Met" = "aliceblue",
                                    "All Assumptions Violated" = "red",
                                    "Low Depth" = "aliceblue",
                                    "Long Timepoint" = "aliceblue"),
                                  c("All Assumptions Met",
                                    "Long Timepoint",
                                    "Low Depth",
                                    "All Assumptions Violated"),
                                  "err_greyout_none.pdf",
                                  FALSE,
                                  c("All Assumptions Violated"))
e4 <- plot_with_assumption_colors(c("All Assumptions Met" = "aliceblue",
                                    "All Assumptions Violated" = "aliceblue",
                                    "Low Depth" = "orange",
                                    "Long Timepoint" = "aliceblue"),
                                  c("All Assumptions Met",
                                    "All Assumptions Violated",
                                    "Long Timepoint",
                                    "Low Depth"),
                                  "err_greyout_low.pdf",
                                  FALSE,
                                  c("Low Depth"))
e5 <- plot_with_assumption_colors(c("All Assumptions Met" = "aliceblue",
                                    "All Assumptions Violated" = "aliceblue",
                                    "Low Depth" = "aliceblue",
                                    "Long Timepoint" = "blue"),
                                  c("All Assumptions Met",
                                    "All Assumptions Violated",
                                    "Low Depth",
                                    "Long Timepoint"),
                                  "err_greyout_long.pdf",
                                  FALSE,
                                  c("Long Timepoint"))


#################
## Other plots ##
#################
ggplot(data = res_all) +
    geom_point(aes(x=depth_reference,
                   y=depth_sample,
                   color=experiment)) +
    geom_vline(xintercept=1307030, color="red", linetype="dashed") +
    geom_hline(yintercept=1307030, color="red", linetype="dashed") +
    geom_abline(slope=1, intercept=0, color="black", linetype="dotted") +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Published Spikeins Have High Variability",
        x = "Depth Reference",
        y = "Depth Sample",
        color="Experiment",
        alpha = "Coverage (X Genome Length)"
    )
ggsave("spikein_consistency.pdf", width = 9, height = 9)

ggplot(data = res_all) +
    geom_point(aes(x=coverage_reference,
                   y=coverage_sample,
                   color=experiment)) +
    geom_vline(xintercept=10, color="red", linetype="dashed") +
    geom_hline(yintercept=10, color="red", linetype="dashed") +
    geom_abline(slope=1, intercept=0, color="black", linetype="dotted") +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    labs(
        title = "Published Spikeins Have High Variability",
        x = "Reference Transcriptome Coverage",
        y = "Sample Transcriptome Coverage",
        color="Experiment",
        alpha = "Coverage (X Transcriptome Length)"
    )
ggsave("coverage_consistency.pdf", width = 9, height = 9)

ggplot(data = res_all) +
    geom_point(aes(x=coverage_sample,
                   ## y=abs(log2(factor_virtual_human)-
                   ##       log2(factor_virtual_drosophila)),
                   y=var_virtual_drosophila,
                   color=experiment), size=3) +
    geom_vline(xintercept=10, color="red", linetype="dotted") +
    ## geom_vline(xintercept=5, color="blue", linetype="dotted") +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16)) +
    ## scale_x_log10() +
    ylim(c(0.8,1.0)) +
    labs(
        title = "Low Sequenced Spike-Ins Have Variable Error Estimates",
        x = "Sample Spike-In Genome Coverage",
        y = "VSI Error Estimate",
        color="Experiment",
        alpha = "Coverage (X Transcriptome Length)"
    )
ggsave("estimate_variability.pdf", width = 9, height = 9)

ggplot(data=res_all) +
    geom_histogram(aes(x=coverage_sample)) +
    ## Assume 30mil deep experiment on human hg38 genome:
    geom_vline(aes(xintercept=10, color="red")) +
    labs(x = 'Average Coverage per Transcriptome Base',
         y = 'Count',
         title = paste0('Spike-In Transcriptome Coverage (',
                        n_expt, " Publications, n=", n_sample, ")")) +
    theme_tufte() +
    theme(text=element_text(family="Helvetica", size=16), legend.position="none")
## print(plotDir)
ggsave("spikein_millionsmapped_coverage.pdf", width=8, height=6)

data.frame(res_all %>% filter(experiment == "Aoi2020nelf"))

######################################################################
### compare_replicates.R ends here
