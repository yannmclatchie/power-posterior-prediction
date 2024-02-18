library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(simstudy)
source("R/normal-location-misspecified/config.R")

# numerical computation of TVD
dnormmix <- function(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                     sigma_3, lambda_1, lambda_2, lambda_3) {
  # compute the density of the gaussian mixture
  return(lambda_1 * dnorm(x, mean = mu_1, sd = sigma_1)
         + lambda_2 * dnorm(x, mean = mu_2, sd = sigma_2)
         + lambda_3 * dnorm(x, mean = mu_3, sd = sigma_3))
}
tvd_integrand <- function(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                          sigma_3, lambda_1, lambda_2, lambda_3,
                          mu_pred, sigma_pred) {
  0.5 * abs(dnorm(x, mean = mu_pred, sd = sigma_pred)
            - dnormmix(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                       sigma_3, lambda_1, lambda_2, lambda_3))
}
tvd_normal_misspecified <- function(mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                                    sigma_3, lambda_1, lambda_2, lambda_3,
                                    mu_pred, sigma_pred) {
  integrate(tvd_integrand, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
            sigma_3, lambda_1, lambda_2, lambda_3,
            mu_pred, sigma_pred, lower = -Inf, upper = Inf)$value
}

# main experiment function
tau_tvd <- function(iter, n, prior, tau, dgp) {
  # simulate data from the true DGP
  dx <- genData(n, dgp$def)
  y <- dx$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # extract mixture variance
  sigma_mix <- dgp$sigma_mix
  
  # compute the posterior parameters
  bar_y <- mean(y)
  sigma_post <- sqrt(1 / ((n * tau) / sigma_mix^2 + 1 / sigma_0^2))
  mu_pred <- sigma_post^2 * (mu_0 / sigma_0^2 + (n * tau * bar_y) / sigma_mix^2)
  sigma_pred <- sqrt(sigma_post^2 + sigma_mix^2)
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_normal_misspecified(mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                                 sigma_3, lambda_1, lambda_2, lambda_3,
                                 mu_pred, sigma_pred) 
  
  # return output
  return(list(tvd = tvd,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# evaluate the elpd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) tau_tvd(iter = iter, 
                                      n = n,
                                      tau = tau,
                                      prior = prior,
                                      dgp = dgp),
       .progress = TRUE) |> # progress bar
  bind_rows()

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# save resutls to csv 
file_name <- paste0("data/normal-location-misspecified-tvd.csv")
write_csv(df, file = file_name)

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
ggsave("./figs/normal-location-misspecified-tvd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_tvd, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-misspecified-tvd.tex")
