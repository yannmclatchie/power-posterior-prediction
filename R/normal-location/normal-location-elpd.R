library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rmutil)
library(reshape2)
library(lemon)
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
  n <- length(y_loo)
  bar_y_i <- mean(y_loo)
  sigma_i <- sqrt(1 / ((n * tau) / sigma_ast^2 + 1 / sigma_0^2))
  mu_i <- sigma_i^2 * (mu_0 / sigma_0^2 + (n * tau * bar_y_i) / sigma_ast^2)
  
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
facet_names <- c(
  `elpd` = "elpd(tau)",
  `tau` = "log2(tau)",
  `2` = "n = 2",
  `10` = "n = 10",
  `100` = "n = 100",
  `1000` = "n = 1000"
)
p_tau_sel <- df |>
  # filter based on prior
  filter(prior == "weak") |>
  dplyr::select(-c("prior")) |>
  # mutate tau to be in log scale
  mutate(tau = log2(tau)) |>
  # choose the best tau by elpd
  group_by(n, iter) |>
  filter(elpd == max(elpd)) |>
  ungroup() |>
  melt(id.vars = c("iter", "n")) |>
  as_tibble() |>
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_rep_grid(n ~ variable, switch = "y", scales = "free_x", 
                 labeller = as_labeller(facet_names), repeat.tick.labels = TRUE) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.line.y=element_blank(),
        axis.line.x=element_line(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
p_tau_sel

# save the plot
ggsave("./figs/normal-location-tau-selection.pdf", width = 5, height = 4)
tex_width <- 5 * 0.8; tex_height = 4 * 0.8
save_tikz_plot(p_tau_sel, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-tau-selection.tex")
