library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/normal-location/config.R")

# numerical computation of TVD between a Gaussian and a Student-t
tvd_integrand <- function(x, mu, sigma, df) {
  0.5 * abs(dnorm(x, mean = mu, sd = sigma)
            - dt(x, df = df))
}
tvd_normal_t <- function(mu, sigma, df) {
  integrate(tvd_integrand, mu, sigma, df,
            lower = -Inf, upper = Inf)$value
}

# main experiment function
tau_tvd <- function(iter, n, prior, tau, df_ast) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # compute the posterior parameters
  bar_y <- mean(y)
  sigma_post <- sqrt(1 / ((n * tau) / 1^2 + 1 / sigma_0^2))
  theta_post <- sigma_post^2 * (mu_0 / sigma_0^2 + (n * tau * bar_y) / 1^2)
  sigma_pred <- sqrt(sigma_post^2 + 1^2)
  
  # compute the plug-in parameters
  sigma_hat <- 1
  theta_hat <- bar_y
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_normal_t(theta_post, sigma_pred, df_ast)
  
  return(list(tvd = tvd,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# load in all datasets
datasets <- readRDS("data/datasets/normal-misspecified.RDS")  

# restrict values of tau
taus <- taus[taus > 0 & is.finite(taus)]

# evaluate the tvd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) tau_tvd(iter = iter, 
                                      n = n,
                                      tau = tau,
                                      prior = prior,
                                      df_ast = df_ast),
       .progress = TRUE) |>
  bind_rows()

# save resutls to csv 
#file_name <- paste0("data/normal-location-tvd.csv")
#write_csv(df, file = file_name)

# read results
#df <- read_csv(file_name)

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# scale the tvd by a function of n
#df <- df |>
#  mutate(tvd = tvd / sqrt(n))

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(tvd_mean = mean(tvd),
            tvd_min = quantile(tvd, probs = 0.05),
            tvd_max = quantile(tvd, probs = 0.95))

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Plot the TVD over iterations
p_tvd_misspecified <- ggplot() +
  geom_line(data = df_100,
            aes(tau, tvd, group = iter), 
            colour = "grey",
            alpha = 0.2) +
  geom_ribbon(data = rdf,
              aes(ymin = tvd_min,
                  ymax = tvd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              linetype = "dotted") +
  geom_line(data = rdf,
            aes(tau, tvd_mean),
            size = 0.75) +
  facet_grid(prior ~ n, scales = "free_y") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  xlab("tau") +
  ylab("TVD") +
  paper_theme
p_tvd_misspecified

# save the plot
ggsave("./figs/normal-location-tvd-misspecified.pdf", width = 5, height = 5 / GR)
my_width <- 1
tex_width <- 5 * my_width; tex_height = 2.5 * my_width
save_tikz_plot(p_tvd_misspecified, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-tvd.tex")
