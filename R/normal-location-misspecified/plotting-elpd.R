library(ggplot2)
library(bayesflow)
library(tidyverse)
library(simstudy)
source("R/normal-location-misspecified/config.R")

#files <- list.files("data/normal-location-misspecified-elpd", full.names = TRUE)
#df <- files %>%
#  map(read_csv) %>% 
#  reduce(rbind)
#write_csv(df, "data/normal-location-misspecified-elpd/normal-location-misspecified-elpd-all.csv")
df <- read_csv("data/normal-location-misspecified-elpd/normal-location-misspecified-elpd-all.csv")

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# Group dataframe by iteration
gdf <- df |>
  group_by(n, tau, prior) |>
  summarize(elpd_mean = mean(elpd))

# Identify optimal tau
best_lines <- gdf |>
  ungroup() |>
  group_by(n, prior) |>
  filter(elpd_mean == max(elpd_mean)) 

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(elpd_min = quantile(elpd, probs = 0.1),
            elpd_max = quantile(elpd, probs = 0.9))

# Plot the elpd over iterations
p_elpd <- ggplot() +
  geom_ribbon(data = rdf,
              aes(ymin = elpd_min,
                  ymax = elpd_max,
                  x = tau),
              colour = "grey",
              alpha = 0.,
              linetype = "dotted") +
  geom_line(data = gdf,
            aes(tau, elpd_mean),
            colour = "black",
            size = 1) +
  geom_vline(xintercept = 1,
             linewidth = 0.25,
             linetype = "dashed") +
  facet_grid(prior ~ n, scales = "free_y") +
  scale_x_continuous(trans = "log10", breaks = 10^seq(-4, 4, length.out = 3)) +
  xlab("tau") +
  ylab("elpd(tau)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none") 
p_elpd

# save the plot
ggsave("./figs/normal-location-misspecified-elpd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_elpd, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-misspecified-elpd.tex")
