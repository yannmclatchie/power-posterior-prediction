library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
source("R/beta-binomial/config.R")

# numerical computation of TVD between the DGP and predictive
tvd_beta_binom <- function(n, p, alpha, beta) {
  # convert alpha and beta to m and s
  m <- alpha / (alpha + beta)
  s <- (alpha + beta)
  
  # compute the numerical summation
  sum(0.5 * abs(dbinom(1:n, n, p)
                - dbetabinom(1:n, n, m, s)))
}

# main experiment function
tau_tvd <- function(iter, n, prior, tau, theta_ast) {
  # simulate data from the true DGP
  y <- rbinom(n, 1, theta_ast)

  # extract prior
  alpha <- prior$a
  beta <- prior$b

  # compute posterior parameters
  x <- sum(y)
  z <- n - x
  alpha_n <- alpha + tau * x
  beta_n <- beta + tau * z
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_beta_binom(n, theta_ast, alpha_n, beta_n)
  
  # compute the LOO-CV elpd
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
                                      theta_ast = theta_ast),
       .progress = TRUE) |>
  bind_rows()

# save resutls to csv 
file_name <- paste0("data/beta-binomial-tvd.csv")
write_csv(df, file = file_name)
#df <- read_csv(file = file_name)

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
  facet_wrap( ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log10", 
                     breaks = 10^seq(-4, 4, length.out = 3)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("tau") +
  ylab("TVD") +
  paper_theme
p_tvd

# save the plot
ggsave("./figs/beta-binomial-tvd.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_tvd, width = tex_width, height = tex_height * 0.75,
               filename = "./tikz/beta-binomial-tvd.tex")
