## Generate supplementary plots for presentations
library("tidyverse")
library("ggthemes")

## Change this for different experiments...
setwd('/home/zach/dowell_lab/virtual_spike_in/dat/')

experiment <- c(rep("Experiment 1" , 2) , rep("Experiment 2" , 2), rep("Experiment 3" , 2))
condition <- rep(c("Control" , "Treatment") , 3)
value <- c(20, 10, 15, 8, 25, 13)
data <- as_tibble(data.frame(experiment,condition,value))

ggplot(data = data, aes(x = experiment, y = value, fill = condition)) +
    geom_bar(position="dodge", stat="identity") +
    labs(x = "Experiment", y = "Read Depth", fill = "Condition") +
    theme_tufte() +
    theme(text = element_text(size=30)) +
    ggsave(paste0(getwd(), "/", "rpkm_example.pdf"), width = 12, height = 9)
