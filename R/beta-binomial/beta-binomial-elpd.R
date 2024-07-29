library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rmutil)
source("R/beta-binomial/config.R")

elpd_loo <- function(iter, n, prior, tau, theta_true) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y

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

# load in all datasets
datasets <- readRDS("data/datasets/beta_binomial.RDS")  

# restrict values of tau
taus <- taus[taus > 0 & is.finite(taus)]

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
#write_csv(df, file = file_name)

# read in previously computed results
df <- read_csv(file = "data/beta-binomial-elpd.csv")

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(elpd_min = quantile(elpd, probs = 0.05),
            elpd_max = quantile(elpd, probs = 0.95),
            elpd_mean = mean(elpd))

# Plot the TVD over iterations
p_elpd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, elpd, group = iter), 
            colour = "grey",
            alpha = 0.15) +
  geom_ribbon(data = rdf,
              aes(ymin = elpd_min,
                  ymax = elpd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              linetype = "dotted") +
  facet_wrap( ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  scale_y_continuous(limits = c(-2, 0)) +
  xlab("tau") +
  ylab("elpd") +
  paper_theme
p_elpd

# save the plot
ggsave("./figs/beta-binomial-elpd.pdf", width = 5, height = 5 / GR)
my_width <- 0.9
tex_width <- 5 * my_width; tex_height = (5 / GR) * my_width
save_tikz_plot(p_elpd, width = tex_width, height = tex_height * 0.75,
               filename = "./tikz/beta-binomial-elpd.tex")
