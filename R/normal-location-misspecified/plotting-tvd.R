library(ggplot2)
library(bayesflow)
library(tidyverse)
library(simstudy)
source("R/normal-location-misspecified/config.R")

#files <- list.files("data/normal-location-misspecified-tvd", full.names = TRUE)
#df <- files %>%
#  map(read_csv) %>% 
#  reduce(rbind)
#write_csv(df, "data/normal-location-misspecified-tvd/normal-location-misspecified-tvd-all.csv")
df <- read_csv("data/normal-location-misspecified-tvd/normal-location-misspecified-tvd-all.csv")

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# Group dataframe by iteration
gdf <- df |>
  group_by(n, tau, prior) |>
  summarize(tvd_mean = mean(tvd))

# Identify optimal tau
best_lines <- gdf |>
  ungroup() |>
  group_by(n, prior) |>
  filter(tvd_mean == max(tvd_mean)) 

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(tvd_min = quantile(tvd, probs = 0.1),
            tvd_max = quantile(tvd, probs = 0.9))

# Plot the TVD over iterations
p_tvd <- ggplot() +
  geom_ribbon(data = rdf,
              aes(ymin = tvd_min,
                  ymax = tvd_max,
                  x = tau),
              colour = "grey",
              alpha = 0.,
              linetype = "dotted") +
  geom_line(data = gdf,
            aes(tau, tvd_mean),
            colour = "black",
            size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_grid(prior ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log10", 
                     breaks = 10^seq(-4, 4, length.out = 3)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("tau") +
  ylab("TVD") +
  paper_theme
p_tvd

# save the plot
ggsave("./figs/normal-location-misspecified-tvd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_tvd, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-misspecified-tvd.tex")
