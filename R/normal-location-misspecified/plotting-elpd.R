library(ggplot2)
library(bayesflow)
library(tidyverse)
library(simstudy)
source("R/normal-location-misspecified/config.R")

## data concatenation
#files <- list.files("data/normal-location-misspecified-elpd", full.names = TRUE)
#df <- files %>%
#  map(read_csv) %>% 
#  reduce(rbind)
#write_csv(df, "data/normal-location-misspecified-elpd/normal-location-misspecified-elpd-all.csv")

# data reading
df <- read_csv("data/normal-location-misspecified-elpd/normal-location-misspecified-elpd-all.csv")

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(elpd_min = quantile(elpd, probs = 0.05),
            elpd_max = quantile(elpd, probs = 0.95))

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Plot the elpd over iterations
p_elpd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, elpd, group = iter), 
            colour = "grey",
            #size = 0.2, 
            alpha = 0.2) +
  geom_ribbon(data = rdf,
              aes(ymin = elpd_min,
                  ymax = elpd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              #size = 0.5,
              linetype = "dotted") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_grid(prior ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  xlab("tau") +
  ylab("LOO-CV elpd") +
  paper_theme
p_elpd

# save the plot
ggsave("./figs/normal-location-misspecified-elpd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_elpd, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-misspecified-elpd.tex")
