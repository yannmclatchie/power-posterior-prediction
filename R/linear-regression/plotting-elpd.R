library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/linear-regression/config.R")

## data concatenation
#files <- list.files("data/linear-regression-elpd", full.names = TRUE)
#df <- files %>%
#  map(read_csv) %>% 
#  reduce(rbind)
#write_csv(df, "data/linear-regression-elpd/linear-regression-elpd-all.csv")

# data reading
df <- read_csv("data/linear-regression/linear-regression-elpd-all.csv")

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(elpd_min = quantile(elpd, probs = 0.05),
            elpd_max = quantile(elpd, probs = 0.95),
            elpd_mean = mean(elpd))

# Restrict the data to only the first hundred realisations
df_100 <- df |> filter(iter <= 50)

# Plot the elpd over iterations
p_elpd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, elpd, group = iter), 
            colour = "grey",
            alpha = 0.15) +
  geom_ribbon(data = rdf,
              aes(ymin = elpd_min,
                  ymax = elpd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              linetype = "dotted") +
  facet_grid(prior ~ n, scales = "free_y") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  xlab("tau") +
  ylab("elpd(tau)") +
  paper_theme
p_elpd

# save the plot
ggsave("./figs/linear-regression-elpd.pdf", width = 5, height = 5 / GR)
my_width <- 0.9
tex_width <- 5 * my_width; tex_height = (5 / GR) * my_width
save_tikz_plot(p_elpd, width = tex_width, height = tex_height,
               filename = "./tikz/linear-regression-elpd.tex")
