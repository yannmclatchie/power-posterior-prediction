library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rmutil)
library(reshape2)
library(lemon)
library(patchwork)
source("R/normal-location/config.R")

# LOO-CV elpd of the normal model
elpd_loo <- function(iter, n, prior, tau, theta_ast, sigma_ast) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # compute the LOO-CV log predictive at each observation
  lpd <- numeric(n)
  for (i in 1:n) {
    lpd[i] <- lpd_loo_i(y, i, mu_0, sigma_0, sigma_ast, tau)
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
  n_loo <- length(y_loo) # this is (n - 1)
  bar_y_i <- mean(y_loo)
  sigma_i <- sqrt(1 / ((n_loo * tau) / sigma_ast^2 + 1 / sigma_0^2))
  mu_i <- sigma_i^2 * (mu_0 / sigma_0^2 + (n_loo * tau * bar_y_i) / sigma_ast^2)
  
  # evaluate the log predictive
  log_pred <- dnorm(x = y_oos, mean = mu_i, 
                    sd = sqrt(sigma_i^2 + sigma_ast^2), 
                    log = TRUE)
  return(log_pred)
}

# load in all datasets
datasets <- readRDS("data/datasets/normal.RDS")

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
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# save resutls to csv 
file_name <- paste0("data/normal-location-elpd.csv")
write_csv(df, file = file_name)

# read results from csv
df <- read_csv(file_name)

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
              colour = "grey",
              alpha = 0.15) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_grid(prior ~ n, scales = "free_y") +
  scale_x_continuous(trans = "log10") +
  xlab("tau") +
  ylab("LOO-CV elpd") +
  paper_theme
p_elpd

# histograms
sel_df <- df |>
  # filter based on prior
  filter(prior == "weak") |>
  dplyr::select(-c("prior")) |>
  # mutate tau to be in log scale
  mutate(tau = log(tau)) |>
  # choose the best tau by elpd
  group_by(n, iter) |>
  filter(elpd == max(elpd)) |>
  ungroup() |>
  melt(id.vars = c("iter", "n")) |>
  as_tibble()
p_tau_sel_small <- sel_df |>
  # filter based on n regime
  filter(n < 100) |>
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_rep_grid(n ~ variable, switch = "y", scales = "free_x", 
                 labeller = as_labeller(facet_names), repeat.tick.labels = FALSE) +
  xlab(NULL) +
  ylab(NULL) +
  theme_sel
p_tau_sel_big <- sel_df |>
  # filter based on n regime
  filter(n >= 100) |> 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_rep_grid(n ~ variable, switch = "y", scales = "free_x", 
                 labeller = as_labeller(facet_names), repeat.tick.labels = FALSE) +
  xlab(NULL) +
  ylab(NULL) +
  theme_sel
# patch plots
p_tau_sel <- p_tau_sel_small + p_tau_sel_big
p_tau_sel

# save the plot
ggsave("./figs/normal-location-tau-selection.pdf", width = 6, height = 6 / GR)
my_width <- 1
tex_width <- 6 * my_width; tex_height = 2.25 * my_width
save_tikz_plot(p_tau_sel, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-tau-selection.tex")
