library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(readr)
source("R/normal-location/config.R")

# restrict tau to reasonable values
taus <- taus[taus > 0 & is.finite(taus)]

## ----
## In expectation

# compute the analytic KL risk
risk <- function(n, tau, sigma_0_2) {
  # compute the analytic risk
  sigma_n_2 <- sigma_0_2 / (1 + n * tau * sigma_0_2)
  risk <- (1 / 2 * log(1 + sigma_n_2)
           + (1 + tau * sigma_n_2) / (2 * (1 + sigma_n_2))
           -1 / 2)
  
  # convert the sigma into a prior name
  if (sigma_0_2 == 1) {prior_name <- "weak"}
  if (sigma_0_2 == 1e+26) {prior_name <- "flat"}
  return(list(risk = risk, n = n, tau = tau, sigma_0_2 = sigma_0_2,
              prior = prior_name))
}

# Evaluate the analytic risk across combinations
combis <- expand.grid(n = ns, tau = taus, sigma_0_2 = sigmas)
avg_df <- combis |>
  pmap(\(n, tau, sigma_0_2) risk(n = n,
                                 tau = tau,
                                 sigma_0_2 = sigma_0_2)) |>
  bind_rows()

# plot the risk over tau
prior_names <- c(
  `1` = "normal(0,1)",
  `1e+26` = "flat",
  `2` = "n = 2",
  `10` = "n = 10",
  `100` = "n = 100"
)
p_risk <- ggplot() +
  geom_line(data = avg_df,
            aes(tau, risk),
            colour = "black",
            size = 1) +
  facet_grid(sigma_0_2 ~ n, scales = "free_y",
             labeller = as_labeller(prior_names)) +
  scale_x_continuous(trans = "log10", breaks = 10^seq(-4, 4, length.out = 3)) +
  xlab("tau") +
  ylab("Risk") +
  paper_theme
p_risk

## ----
## Single dataset

kld_tau <- function(iter, n, prior, tau, theta_ast, sigma_ast) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # compute the posterior parameters
  bar_y <- mean(y)
  sigma_post <- sqrt(1 / ((n * tau) / sigma_ast^2 + 1 / sigma_0^2))
  theta_post <- sigma_post^2 * (mu_0 / sigma_0^2 + (n * tau * bar_y) / sigma_ast^2)
  sigma_pred <- sqrt(sigma_post^2 + sigma_ast^2)
  
  # compute the KL
  kld <- (log(sigma_pred / sigma_ast) 
          + (sigma_ast^2 + (theta_ast - theta_post)^2) / (2 * sigma_pred^2) 
          - 1/2)
  
  # return results
  return(list(kld = kld,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# load in all datasets
datasets <- readRDS("data/datasets/normal.RDS")  

# evaluate the tvd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) kld_tau(iter = iter, 
                                      n = n,
                                      tau = tau,
                                      prior = prior,
                                      theta_ast = theta_ast,
                                      sigma_ast = sigma_ast),
       .progress = TRUE) |>
  bind_rows()

# fix ordering of priors
avg_df$prior = factor(avg_df$prior, levels=c('weak', 'flat'))
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# scale by sqrt n
df <- df |>
  mutate(kld = n * kld)
avg_df <- avg_df |>
  mutate(risk = n * risk)

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(kld_min = quantile(kld, probs = 0.05),
            kld_max = quantile(kld, probs = 0.95))

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Plot the TVD over iterations
p_kld <- ggplot() +
  geom_line(data = df_100,
            aes(tau, kld, group = iter), 
            colour = "grey",
            alpha = 0.2) +
  geom_ribbon(data = rdf,
              aes(ymin = kld_min,
                  ymax = kld_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              linetype = "dotted") +
  geom_line(data = avg_df,
            aes(tau, risk),
            colour = "black",
            size = 0.75) +
  facet_grid(prior ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  scale_y_continuous(limits = c(0, 4)) +
  xlab("tau") +
  ylab("n KLD") +
  paper_theme
p_kld

# save the plot
ggsave("./figs/normal-location-kld.pdf", width = 5, height = 5 / GR)
my_width <- 1
tex_width <- 5 * my_width; tex_height = 2.5 * my_width
save_tikz_plot(p_kld, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-kld.tex")
