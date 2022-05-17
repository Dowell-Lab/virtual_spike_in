### gen_comparison_plots.r --- Generate Comparison Plots
##
## Filename: gen_comparison_plots.r
## Author: Zach Maas
## Created: Fri Apr  8 10:02:17 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## This file contains code to generate plots comparing the results of
## Dukler2017nascent using both the endogenous spike-in and the
## computational virtual spike-in method.
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
library('forcats')
library('ggplot2')
library('ggthemes')
library('plyr')
library('extrafont')
loadfonts(device="all")

## Root Analysis Directory
baseDir <- '/Users/zachmaas/Dropbox/phd/research/dna_lab/virtual_spike_in'

## Load Normalization Factors
factors_spikein <- read_delim(
    paste0(baseDir, '/dukler_out/spikein/virtual_size_factors.txt'),

    col_names = c('comparison', 'norm_factor_spikein', 'norm_var_spikein')) %>%
    mutate(comparison = str_remove_all(comparison, '.sorted')) %>%
    separate(comparison, c("reference", "sample"), "_vs_", FALSE)
factors_virtual <- read_delim(
    paste0(baseDir, '/dukler_out/virtual/virtual_size_factors.txt'),
    col_names = c('comparison', 'norm_factor_virtual', 'norm_var_virtual')) %>%
    mutate(comparison = str_remove_all(comparison, '.sorted')) %>%
    separate(comparison, c("reference", "sample"), "_vs_", FALSE)
factors_merged <- full_join(x=factors_spikein, y=factors_virtual, by='comparison') %>%
    separate(comparison, c("reference", "sample"), "_vs_", FALSE)


## Load Count Data
counts_spikein <- read_delim(
    paste0(baseDir, '/dukler_out/spikein/counts/counts_merged.txt'))
counts_virtual <- read_delim(
    paste0(baseDir, '/dukler_out/virtual/counts/counts_merged.txt'))
colnames(counts_spikein) <- gsub('.sorted.sorted.bam',
                                 '', colnames(counts_spikein), fixed=TRUE)
colnames(counts_virtual) <- gsub('.sorted.sorted.bam',
                                 '', colnames(counts_spikein), fixed=TRUE)

## Run Regressions (Endogenous)
## TODO Make sure the A/B is correct here
coefs <- c()
for (row in 1:nrow(factors_virtual)) {
    reference <- factors_virtual[row, "reference"]
    sample <- factors_virtual[row, "sample"]
    formula <- paste0(reference, ' ~ ', sample, collapse=" ")
    reg <- lm(formula, data = counts_virtual)
    coef <- reg$coefficients[2]
    coefs <- c(coefs, coef)
}
factors_virtual$lm_factor_virtual <- coefs

## Run Regressions (Exogenous)
coefs <- c()
for (row in 1:nrow(factors_spikein)) {
    reference <- factors_spikein[row, "reference"]
    sample <- factors_spikein[row, "sample"]
    formula <- paste0(reference, ' ~ ', sample, collapse=" ")
    print(formula)
    reg <- lm(formula, data = counts_spikein)
    coef <- reg$coefficients[2]
    coefs <- c(coefs, coef)
}
factors_virtual$lm_factor_spikein <- coefs

## SRR/Timepoint Lists
srrs <- c("SRR5364303","SRR5364304",
          "SRR5364305","SRR5364306",
          "SRR5364307","SRR5364308",
          "SRR5364309","SRR5364310",
          "SRR5364311","SRR5364312",
          "SRR5364313","SRR5364314")
timepoints <- c("00min_Rep1","00min_Rep2",
                "10min_Rep1","10min_Rep2",
                "20min_Rep1","20min_Rep2",
                "40min_Rep1","40min_Rep2",
                "60min_Rep1","60min_Rep2",
                "160min_Rep1","160min_Rep2")
## Original version...
reads_mapped_c <- c(31057828, 32284626,
                    30856790, 28626647,
                  28646728, 25865333,
                  31928691, 35265731,
                  36406958, 27779495,
                  34118092, 36305470)
## Fixed version from hisat2
reads_mapped <- c(41341 + 188966, 48240 + 171766,
                  50000 + 117447, 52709 + 129299,
                   52608 + 108062, 43874 + 110758,
                  38029 + 101125, 61073 + 104243,
                  58209 + 92171, 41553 + 74251,
                  47768 + 93566, 61827 + 111719)
## Switch labels to timepoints
factors_virtual$reference_mapped <- as.numeric(
    mapvalues(factors_virtual$reference,
              from=srrs,
              to=reads_mapped))
factors_virtual$reference <- mapvalues(factors_virtual$reference,
                                       from=srrs,
                                       to=timepoints)
factors_virtual$sample_mapped <- as.numeric(
    mapvalues(factors_virtual$sample,
              from=srrs,
              to=reads_mapped))
factors_virtual$sample <- mapvalues(factors_virtual$sample,
                                    from=srrs,
                                    to=timepoints)
factors_virtual <- factors_virtual %>%
    mutate(comparison=fct_rev(fct_inorder(paste0(str_replace(reference, "_", " "),
                                                 " vs ",
                                                 str_replace(sample, "_", " ")))),
           sample=fct_rev(fct_inorder(str_replace(sample, "_", " "))),
           mapped_ratio = reference_mapped / sample_mapped)

## TODO Bump up to 60 minute data set (180kb)
## TODO Show region as onset of apoptosis
## Generate plots
plotDir <- paste0(baseDir, '/dukler_out/plots/')
for (row in 1:nrow(factors_virtual)) {
    comparison <- factors_virtual[row, "comparison"]
    sample <- factors_virtual[row, "sample"]
    ## Grab data
    norm_factor <- log2(as.numeric(factors_virtual[row, "norm_factor_virtual"]))
    norm_var <- log2(as.numeric(factors_virtual[row, "norm_var_virtual"]))
    lm_factor_mahat <- log2(as.numeric(factors_virtual[row, "lm_factor_virtual"]))
    mapped_ratio <- log2(as.numeric(factors_virtual[row, "mapped_ratio"]))
    ## Calculate deviation
    virtual_deviation <- round(abs(mapped_ratio - norm_factor), digits=3)
    mahat_deviation <- round(abs(mapped_ratio - lm_factor_mahat), digits=3)
    virtual_label <- paste0("VSI Distribution, difference = ", virtual_deviation)
    mahat_label <- paste0("Linear Regression, difference = ", mahat_deviation)
    ## Make plot
    ggplot() +
        stat_function(fun = dnorm, n = 200,
                      args = list(mean = norm_factor,
                                  sd = sqrt(norm_var)),
                      aes(color="black")) +
        geom_vline(aes(xintercept = lm_factor_mahat,
                       color="blue")) +
        geom_vline(aes(xintercept = mapped_ratio,
                       color="green")) +
        geom_vline(xintercept = 0, linetype="dotted", color="red") +
        xlim(c(-2, 2)) +
        scale_color_manual(name="Method",
                           values = c("black" = "black",
                                      "blue" = "blue",
                                      "green" = "green"),
                           labels=c(virtual_label,
                                    mahat_label,
                                    "Mapped Ratio (Exogenous)")) +
        labs(title = comparison,
             x = "log2(Normalization Factor)",
             y = "Relative Proportion") +
        theme_tufte() +
        theme(text=element_text(family="Helvetica"))
    ## scale_fill_discrete(name = "Dose", labels = c("A", "B", "C")) +
    ggsave(paste0(plotDir, sample, "_comparison.pdf"),
           width=8, height=6)
}

ggplot(data=factors_virtual) +
    geom_pointrange(aes(x=sample, y=log2(norm_factor_virtual),
                        ymin=log2(norm_factor_virtual)-log2(norm_var_virtual),
                        ymax=log2(norm_factor_virtual)+log2(norm_var_virtual),
                        color="black")) +
    geom_pointrange(aes(x=sample, y=log2(lm_factor_virtual),
                        ymin=log2(lm_factor_virtual),
                        ymax=log2(lm_factor_virtual),
                        color="blue")) +
    geom_pointrange(aes(x=sample, y=log2(mapped_ratio),
                        ymin=log2(mapped_ratio),
                        ymax=log2(mapped_ratio),
                        color="green")) +
    geom_hline(yintercept = 0, linetype="dotted", color="red") +
    scale_color_manual(name="Method",
                       values = c("black" = "black",
                                  "blue" = "blue",
                                  "green" = "green"),
                       labels=c("VSI",
                                "Linear Regression",
                                "Mapped Ratio (Exogenous)")) +
    ylim(c(-2.5, 2.5)) +
    labs(title = "Comparisons to 00min Rep1",
         y = "log2(Normalization Factor)",
         x = "Sample Comparison") +
    theme_tufte() +
    theme(legend.position="bottom", legend.justification="right",
          text=element_text(family="Helvetica")) +
    coord_flip()
ggsave(paste0(plotDir, "dukler_comparison_all.pdf"),
       width=6, height=8)

######################################################################
### gen_comparison_plots.r ends here
