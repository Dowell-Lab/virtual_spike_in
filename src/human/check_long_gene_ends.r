## Script to run DESeq2
## Maintainer: Zachary Maas <zama8258@colorado.edu>

#######################
## Preliminary Setup ##
#######################

library("tidyverse")
library("ggthemes")
library("pbapply")
library("gridExtra")
library("DESeq2")

## Change this for different experiments...
setwd('/home/zach/dowell_lab/virtual_spike_in/dat/')
num_reps <- 2

## To filter out only the genes we've used and convert to common ID
idConvert <- read_delim("refseq_to_common_id.txt",
                        col_names = c("rowname", "common"), delim = " ")
## Experimental condition
condition <- factor(c(rep("treated", num_reps), rep("untreated", num_reps)), levels=c("untreated", "treated"))

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
full_seqdata <- read.delim("counts_fix.txt", stringsAsFactors = FALSE, header = TRUE,
                           row.names = 1)
## Then, we filter to only include the resolved isoforms
full_df <- as_tibble(rownames_to_column(full_seqdata))

## For DROSOPHILA
if (!is.null(full_df$rowname)) {
    full_df$common <- full_df$rowname
    full_df$rowname <- NULL
}

full_dt <- data.frame(full_df)
rownames(full_dt) <- full_dt$common
full_dt$common <- NULL
full_dt <- as_tibble(full_dt)

full_dt <- full_dt %>%
    mutate(wt_mean = (rowMeans(full_dt[c('wt_U2OS_b1.sorted.bam',
                                         'wt_U2OS_b2.sorted.bam')])),
           wt1h_mean = (rowMeans(full_dt[c('wt_1hAs_U2OS_b1.sorted.bam',
                                           'wt_1hAs_U2OS_b2.sorted.bam')])),
           dg_mean = (rowMeans(full_dt[c('dG3BP_U2OS_b1.sorted.bam',
                                         'dG3BP_U2OS_b2.sorted.bam')])),
           dg1h_mean = (rowMeans(full_dt[c('dG3BP_1hAs_U2OS_b1.sorted.bam',
                                           'dG3BP_1hAs_U2OS_b2.sorted.bam')]))) %>%
    subset(select = c('wt_mean', 'wt1h_mean', 'dg_mean', 'dg1h_mean'))


## Use linear regression to check for baseline skew between samples
library("speedglm")
reg <- function(df) {
    slopes <- c(
        summary(lm(dg_mean   ~ 0 + wt_mean, data=df  ))$coefficients[1],
        summary(lm(wt_mean   ~ 0 + wt1h_mean, data=df))$coefficients[1]
    )
    return(slopes)
}

## NOTE: This step is not required, it's just for checking the quality
## of results.

## Perform monte carlo subsampling to verify those results. Zhang, P.
## (1993). Model Selection Via Muiltfold Cross Validation. Ann. says
## that using N^2 will give us close to optimal results. That will
## take too long on my laptop, so we'll just do 100k trials, which is
## about 10% of that value.
sample_slopes <- pbreplicate(100000, reg(sample_frac(full_dt, 0.15, replace = TRUE)))
slope_dist <- as_tibble(t(sample_slopes))
colnames(slope_dist) <- c("dg_v_wt",
                          "wt_v_wt1h")

## Mapped data stats
## mapstats <- tibble(dg1h_1 = 906583, dg1h_2 = 821856, dg_1 = 1299042,
##                    dg_2 = 821856, wt1h_1 = 1784247, wt1h_2 = 334378,
##                    wt_1 = 1988813, wt_2 = 2376890) %>%
mapstats <- tibble(dg1h_1 = 56958302, dg1h_2 = 69950273, dg_1 = 58282107,
                   dg_2 = 71052598, wt1h_1 = 91587100, wt1h_2 = 66604306,
                   wt_1 = 85780194, wt_2 = 93689147) %>%
    mutate(dg1h = mean(dg1h_1, dg1h_2),
           dg = mean(dg_1, dg_2),
           wt1h = mean(wt1h_1, wt1h_2),
           wt = mean(wt_1, wt_2)) %>%
    mutate(dg_v_dg1h = dg / dg1h,
           dg_v_wt = dg / wt,
           dg_vs_wt1h = dg / wt1h,
           dg1h_v_wt = dg1h / wt,
           dg1h_v_wt1h = dg1h / wt1h,
           wt_v_wt1h = wt / wt1h) %>%
    subset(select = c(dg_v_dg1h, dg_v_wt, dg_vs_wt1h,
                      dg1h_v_wt, dg1h_v_wt1h, wt_v_wt1h))


## Generate some distribution plots for the cross validation
dg_v_wt_pl <- ggplot() +
    geom_density(data = slope_dist, aes(x = dg_v_wt)) +
    geom_vline(data = mapstats, aes(color = "red", xintercept = dg_v_wt), show.legend=FALSE) +
    labs(title = "Mutant at 0hr vs Wild Type at 0hr", x = "Normalization Factor",
         y = "Relative Proportion") + theme_tufte() + xlim(0, 2)
wt_v_wt1h_pl <- ggplot() +
    geom_density(data = slope_dist, aes(x = wt_v_wt1h)) +
    geom_vline(data = mapstats, aes(color = "red", xintercept = wt_v_wt1h), show.legend=FALSE) +
    labs(title = "Wild Type at 0hr vs Wild Type with 1hr Arsenic", x = "Normalization Factor",
         y = "Relative Proportion") + theme_tufte() + xlim(0, 2)
arr <- grid.arrange(dg_v_wt_pl, wt_v_wt1h_pl,
                    nrow = 1,
                    top = "Long-Gene End Normalization Distribution vs Spike-In Normalization",
                    bottom = "Vertical Lines use Spike-In, Distribution uses Ends of Long Genes")
ggsave("virtual_spike_in_comparison.pdf", arr, width = 10, height = 5, device = "pdf")
