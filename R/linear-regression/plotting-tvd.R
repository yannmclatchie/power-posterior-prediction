library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/linear-regression/config.R")

## data concatenation
#files <- list.files("data/linear-regression-tvd", full.names = TRUE)
#df <- files %>%
#  map(read_csv) %>% 
#  reduce(rbind)
#write_csv(df, "data/linear-regression-tvd/linear-regression-tvd-all.csv")

# data reading
df <- read_csv("data/linear-regression/linear-regression-tvd-all.csv")

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(tvd_min = quantile(tvd, probs = 0.05),
            tvd_max = quantile(tvd, probs = 0.95),
            tvd_mean = mean(tvd))

# Choose 50 random iterations
df_100 <- df |> filter(iter %in% sample(unique(df$iter), 50))

# Plot the tvd over iterations
p_tvd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, tvd, group = iter), 
            colour = "grey",
            alpha = 0.15) +
  geom_ribbon(data = rdf,
              aes(ymin = tvd_min,
                  ymax = tvd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              linetype = "dotted") +
  facet_grid(prior ~ n, scales = "free") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  xlab("tau") +
  ylab("tvd") +
  paper_theme
p_tvd

# save the plot
ggsave("./figs/linear-regression-tvd.pdf", width = 5, height = 5 / GR)
my_width <- 0.9
tex_width <- 5 * my_width; tex_height = (5 / GR) * my_width
save_tikz_plot(p_tvd, width = tex_width, height = tex_height,
               filename = "./tikz/linear-regression-tvd.tex")

# use just one box for the main body
df_single <- df |>
  filter(n == 50, prior == "weak")
df_single_100 <- df_single |> filter(iter %in% sample(unique(df$iter), 50))
rdf_single <- df_single |>
  group_by(tau, n, prior) |>
  summarize(tvd_min = quantile(tvd, probs = 0.05),
            tvd_max = quantile(tvd, probs = 0.95),
            tvd_mean = mean(tvd))

# make single-facet plot
p_tvd_single <- ggplot() +
  geom_line(data = df_single_100,
            aes(tau, tvd, group = iter), 
            colour = "grey",
            #size = 0.2, 
            alpha = 0.15) +
  geom_ribbon(data = rdf_single,
              aes(ymin = tvd_min,
                  ymax = tvd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              #size = 0.5,
              linetype = "dotted") +
  facet_wrap(~ n, scales = "free") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  xlab("tau") +
  ylab("tvd") +
  paper_theme
p_tvd_single

# save the plot
ggsave("./figs/linear-regression-tvd.pdf", width = 5, height = 5 / GR)
my_width <- 0.8
scaling <- 0.75
tex_width <- 5 * my_width * 0.6; tex_height = (5 / GR) * my_width
tex_width <- tex_width * scaling; tex_height = tex_height * scaling
save_tikz_plot(p_tvd_single, width = tex_width, height = tex_height,
               filename = "./tikz/linear-regression-tvd-single.tex")
