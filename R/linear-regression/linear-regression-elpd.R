library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/linear-regression/config.R")

# main experiment function
elpd_loo <- function(iter, n, prior, tau, theta_ast, sigma_ast) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  X <- data$X
  y <- data$y
  
  # extract prior
  mu_0 <- prior$mu_0
  Sigma_0 <- prior$Sigma_0
  
  # compute the LOO-CV log predictive at each observation
  lpd <- numeric(n)
  for (i in 1:n) {
    lpd[i] <- lpd_loo_i(y, X, i, Sigma_0, sigma_ast, tau)
  }
  
  # compute the TVD between the posterior predictive and the true DGP
  return(list(elpd = mean(lpd),
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}
lpd_loo_i <- function(y, X, i, Sigma_0, sigma_ast, tau) {
  # produce the deleted data
  y_loo <- y[-i]; X_loo <- as.matrix(X[-i,])
  y_oos <- y[i]; X_oos <- as.vector(X[i,])
  
  # fix dimensions in the n = 2 case
  if (length(y) == 2) {
    X_loo <- t(X_loo)
  }
  
  # compute the posterior parameters
  Sigma_n <- tau / sigma_ast^2 * t(X_loo) %*% X_loo + solve(Sigma_0)
  theta_n <- 1 / sigma_ast^2 * Sigma_n %*% t(X_loo) %*% y_loo

  # compute the posterior predictive parameters
  mu_pred <- t(theta_n) %*% X_oos
  sigma_pred <- sqrt(t(X_oos) %*% Sigma_n %*% X_oos + sigma_ast^2)
  print(c(n, tau, mu_pred, sigma_pred))
  
  # evaluate the leave-one-out predictive
  lpd <- dnorm(x = y_oos, mean = mu_pred, sd = sigma_pred, log = TRUE)
  return(lpd)
}

# load in all datasets
datasets <- readRDS("data/datasets/linear_regression.RDS")

# Evaluate the elpd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) elpd_loo(iter = iter, 
                                       n = n,
                                       tau = tau,
                                       prior = prior,
                                       theta_ast = theta_ast,
                                       sigma_ast = sigma_ast),
       .progress = TRUE) |>
  bind_rows()

# fix ordering of priors
#df$prior = factor(df$prior, levels=c('weak', 'flat'))

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(elpd_min = quantile(elpd, probs = 0.05),
            elpd_max = quantile(elpd, probs = 0.95))

# Restrict the data to only the first hundred realisations
df_100 <- df |> filter(iter <= 100)

# Plot the elpd over iterations
p_elpd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, elpd, group = iter), 
            colour = "grey",
            #size = 0.2, 
            alpha = 0.15) +
  geom_ribbon(data = rdf,
              aes(ymin = elpd_min,
                  ymax = elpd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              #size = 0.5,
              linetype = "dotted") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~ n, scales = "free") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  #scale_y_continuous(trans = "pseudo_log") +
  xlab("tau") +
  ylab("elpd(tau)") +
  paper_theme
p_elpd

