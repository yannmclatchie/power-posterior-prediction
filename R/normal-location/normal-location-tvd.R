library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/normal-location/config.R")

# numerical computation of TVD between two normals
tvd_integrand <- function(x, mu_1, sigma_1, mu_2, sigma_2) {
  0.5 * abs(dnorm(x, mean = mu_1, sd = sigma_1)
            - dnorm(x, mean = mu_2, sd = sigma_2))
}
tvd_normal <- function(mu_1, sigma_1, mu_2, sigma_2) {
  integrate(tvd_integrand, mu_1, sigma_1, mu_2, sigma_2,
            lower = -Inf, upper = Inf)$value # rel.tol = 1e-9
}

# main experiment function
tau_tvd <- function(iter, n, prior, tau, theta_ast, sigma_ast) {
  # simulate data from the true DGP
  y <- rnorm(n, theta_ast, sigma_ast)
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # compute the posterior parameters
  bar_y <- mean(y)
  sigma_post <- sqrt(1 / ((n * tau) / sigma_ast^2 + 1 / sigma_0^2))
  theta_post <- sigma_post^2 * (mu_0 / sigma_0^2 + (n * tau * bar_y) / sigma_ast^2)
  sigma_pred <- sqrt(sigma_post^2 + sigma_ast^2)
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_normal(theta_ast, sigma_ast, theta_post, sigma_pred)
    return(list(tvd = tvd,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# evaluate the tvd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) tau_tvd(iter = iter, 
                                      n = n,
                                      tau = tau,
                                      prior = prior,
                                      theta_ast = theta_ast,
                                      sigma_ast = sigma_ast),
       .progress = TRUE) |>
  bind_rows()

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# save resutls to csv 
file_name <- paste0("data/normal-location-tvd.csv")
write_csv(df, file = file_name)

df <- read_csv(file_name)

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
  xlab("tau") +
  ylab("TVD") +
  paper_theme
p_tvd

# save the plot
ggsave("./figs/normal-location-tvd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_tvd, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-tvd.tex")

# compute the Hellinger bound
tau_hel_sq <- function(n, tau, k0, theta_ast, sigma_ast) {
  (1 - sqrt((4 * sqrt(1 + 1 / (n * tau + k0))) 
           / (2 * (2 + 1 / (n * tau + k0)) + n * tau^2 / (n * tau + k0)^2))
   * exp(-(k0^2 * theta_ast^2) / 
           (2 * sigma_ast^2 * (2 * (2 * (n * tau + k0) + 1) 
                               * (n * tau + k0) + n * tau^2))))
}
hel_bounds <- function(n, tau, prior, theta_ast, sigma_ast) {
  # compute the k0 factor
  k0 = sigma_ast^2 / prior$sigma^2

  # return the bounds
  return(list(n = n,
              tau = tau,
              prior = prior$name,
              lower = tau_hel_sq(n, tau, k0, theta_ast, sigma_ast),
              upper = sqrt(2) * sqrt(tau_hel_sq(n, tau, k0, theta_ast, sigma_ast))))
}
combis <- expand.grid(n = ns, tau = taus, prior = priors)
hel_df <- combis |>
  pmap(\(n, tau, prior) hel_bounds(n = n,
                                   tau = tau,
                                   prior = prior,
                                   theta_ast = theta_ast,
                                   sigma_ast = sigma_ast)) |>
  bind_rows()
  
# add Hellinger bounds to plot
p_tvd +
  geom_ribbon(data = hel_df,
              aes(tau, ymax = upper, ymin = lower),
              fill = NA,
              colour = "red")

