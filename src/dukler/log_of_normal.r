## Small experiment for the log of normal distributions

library('tidyverse')
library('MASS')

dist <- rnorm(10000,5,0.5)
log_dist <- log2(dist)

fitdistr(log_dist, "normal")

ggplot() +
    geom_density(aes(x=log_dist))
