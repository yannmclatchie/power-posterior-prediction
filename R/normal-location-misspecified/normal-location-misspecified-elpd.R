library(simstudy)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(bayesflow)
source("R/normal-location-misspecified/config.R")

# compute the LOO-CV elpd 
elpd_loo <- function(iter, n, prior, tau, dgp) {
  # simulate data from the true DGP
  dx <- genData(n, dgp$def)
  y <- dx$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # extract mixutre variance
  sigma_mix <- dgp$sigma_mix
  
  # compute the LOO-CV log predictive at each observation
  lpd <- numeric(n)
  for (i in 1:n) {
    lpd[i] <- lpd_loo_i(y, i, mu_0, sigma_0, sigma_mix, tau)
  }
  
  # compute the LOO-CV elpd
  return(list(elpd = sum(lpd) / n,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}
lpd_loo_i <- function(y, i, mu_0, sigma_0, sigma_ast, tau) {
  # produce the deleted data
  y_loo <- y[-i]
  y_oos <- y[i]
  
  # compute posterior parameters
  n <- length(y)
  bar_y_i <- mean(y_loo)
  sigma_i <- sqrt(1 / (((n - 1) * tau) / sigma_ast^2 + 1 / sigma_0^2))
  mu_i <- sigma_i^2 * (mu_0 / sigma_0^2 + ((n - 1) * tau * bar_y_i) / sigma_ast^2)
  
  # evaluate the log predictive
  log_pred <- dnorm(x = y_oos, mean = mu_i, 
                    sd = sqrt(sigma_i^2 + sigma_ast^2), log = TRUE)
  return(log_pred)
}

# Evaluate the elpd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) elpd_loo(iter = iter, 
                                       n = n,
                                       tau = tau,
                                       prior = prior,
                                       dgp = dgp),
       .progress = TRUE) |>
  bind_rows()

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# save resutls to csv 
file_name <- paste0("data/normal-location-misspecified-elpd.csv")
write_csv(df, file = file_name)

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
