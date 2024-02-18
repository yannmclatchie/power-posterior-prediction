library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rmutil)
source("R/beta-binomial/config.R")

elpd_loo <- function(iter, n, prior, tau, theta_true) {
  # simulate data from the true DGP
  y <- rbinom(n, 1, theta_true)
  
  # extract prior
  alpha <- prior$a
  beta <- prior$b
  
  # compute the LOO-CV log predictive at each observation
  lpd <- numeric(n)
  for (i in 1:n) {
    lpd[i] <- lpd_loo_i(y, i, alpha, beta, tau)
  }
  
  # compute the LOO-CV elpd
  return(list(elpd = sum(lpd) / n,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

lpd_loo_i <- function(y, i, alpha, beta, tau) {
  # produce the deleted data
  y_loo <- y[-i]
  y_oos <- y[i]
  
  # compute posterior parameters
  n <- length(y_loo)
  x <- sum(y_loo)
  z <- n - x
  alpha_n <- alpha + tau * x
  beta_n <- beta + tau * z
  
  # evaluate the log predictive
  m <- alpha_n / (alpha_n + beta_n)
  s <- (alpha_n + beta_n)
  log_pred <- dbetabinom(y_oos, 1, m, s, log = TRUE)
  return(log_pred)
}

# Evaluate the elpd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) elpd_loo(iter = iter, 
                                       n = n,
                                       tau = tau,
                                       prior = prior,
                                       theta_true = theta_true),
       .progress = TRUE) |>
  bind_rows()

# save resutls to csv 
file_name <- paste0("data/beta-binomial-elpd.csv")
write_csv(df, file = file_name)

# Group dataframe by iteration
gdf <- df |>
  group_by(tau, n, prior) |>
  summarize(elpd_mean = mean(elpd))

# Plot the mean results
gdf |> 
  ggplot(aes(x = tau, y = elpd_mean)) + 
  geom_line() +
  facet_wrap(~n, scales = "fixed")

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
  facet_wrap(~n, scales = "fixed", nrow = 1) +
  scale_x_continuous(trans = "log10", breaks = 10^seq(-4, 4, length.out = 3)) +
  coord_cartesian(ylim = c(-3, 0)) +
  scale_y_continuous(limits = c()) +
  xlab("tau") +
  ylab("elpd(t)") +
  paper_theme
p_elpd

# save the plot
ggsave("./figs/beta-binomial-elpd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_elpd, width = tex_width, height = tex_height * 0.75,
               filename = "./tikz/beta-binomial-elpd.tex")

