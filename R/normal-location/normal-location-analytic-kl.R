library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggmagnify)
library(directlabels)
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
  if (sigma_0_2 == sqrt(0.2)) {prior_name <- "informative"}
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

## ---
## Pre-asymptotic

# define pre-asymptotic regime experiments
n_range <- seq(1, 500)
taus <- c(1, 0.5, 0.05)

# evaluate the analytic risk across combinations
asymp_combis <- expand.grid(n = n_range, tau = taus)
asymp_df <- asymp_combis |>
  pmap(\(n, tau) risk(n = n,
                      tau = tau,
                      sigma_0_2 = sqrt(0.2))) |>
  bind_rows()

# plot the pre-asymptotic regime
asymp_p <- asymp_df |> 
  mutate(risk = n * risk,
         tau = as.factor(tau)) |>
  ggplot(aes(n, risk, colour = tau, linetype = tau)) +
  geom_line() +
  scale_colour_manual(values = c("black", "grey", "black")) +
  scale_linetype_manual(values = c("solid", "solid", "dotdash")) +
  #scale_colour_manual(values = c("#4477AA", "grey", "#996633")) +
  paper_theme
  #theme_bw() +
  #theme(legend.position = "bottom")

# add magnification
from <- c(xmin = 0, xmax = 20, ymin = 0, ymax = 0.5)
to <- c(xmin = 200, xmax = 450, ymin = 0, ymax = 0.4)
asymp_p <- asymp_p + 
  geom_magnify(from = from, to = to, axes = "xy")

# add direct labels
asymp_p <- direct.label(asymp_p,"maxvar.points")
asymp_p

# save the plot
my_width <- 1
tex_width <- 5 * my_width; tex_height = 2.5 * my_width
save_tikz_plot(asymp_p, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-pre-asymp.tex")

# define alternative pre-asymptotic regime experiments
n_range <- c(seq(1, 50), seq(100, 1e4))
taus <- c(1, 0.5, 0.1)

# evaluate the analytic risk across combinations
asymp_combis <- expand.grid(n = n_range, tau = taus)
asymp_df <- asymp_combis |>
  pmap(\(n, tau) risk(n = n,
                      tau = tau,
                      sigma_0_2 = sqrt(0.2))) |>
  bind_rows()

# alternative plot
alt_asymp_p <- asymp_df |> 
  mutate(risk = n * risk,
         tau = as.factor(tau),
         data_regime = n >= 50) |>
  ggplot(aes(n, risk, colour = tau, linetype = tau)) +
  geom_line() +
  facet_wrap(~data_regime, scales = "free", ncol=2) +
  scale_colour_manual(values = c("black", "grey", "black")) +
  scale_linetype_manual(values = c("solid", "solid", "dotdash")) +
  paper_theme + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
alt_asymp_p

# save the plot
my_width <- 1
tex_width <- 5 * my_width; tex_height = 2.5 * my_width
save_tikz_plot(alt_asymp_p, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-pre-asymp-alt.tex")

## ---
## temperature taken at rate 1/n

# compute the analytic KL risk in terms of an alpha
risk_alpha <- function(n, alpha, sigma_0_2) {
  # compute tau
  tau <- alpha / (alpha + n)
  
  # compute the analytic risk
  sigma_n_2 <- sigma_0_2 / (1 + n * tau * sigma_0_2)
  risk <- (1 / 2 * log(1 + sigma_n_2)
           + (1 + tau * sigma_n_2) / (2 * (1 + sigma_n_2))
           -1 / 2)
  
  # convert the sigma into a prior name
  if (sigma_0_2 == sqrt(0.2)) {prior_name <- "informative"}
  if (sigma_0_2 == 1) {prior_name <- "weak"}
  if (sigma_0_2 == 1e+26) {prior_name <- "flat"}
  return(list(risk = risk, n = n, alpha = alpha, tau = tau, 
              sigma_0_2 = sigma_0_2, prior = prior_name))
}

# define ranges for the experiment
n_range <- seq(1, 500)
alphas <- rexp(n = 6, rate = 1)

# evaluate the analytic risk across combinations
coarse_combis <- expand.grid(n = n_range, alpha = alphas)
coarse_df <- coarse_combis |>
  pmap(\(n, alpha) risk_alpha(n = n,
                              alpha = alpha,
                              sigma_0_2 = sqrt(0.2))) |>
  bind_rows()

# plot the (non-)convergence of the coarsened posterior
coarse_p <- coarse_df |> 
  mutate(risk = risk,
         alpha = paste0("alpha = ", as.factor(round(alpha, digits = 3)))) |>
  ggplot(aes(n, risk, linetype = alpha, fill = alpha)) +
  geom_line() +
  geom_hline(yintercept = 0, size = 0.3) + 
  scale_colour_manual(values = rep("black", times = 6)) +
  paper_theme +
  ylim(0, max(coarse_df$risk)) + 
  xlim(0, 550)
coarse_p <- direct.label(coarse_p,"last.qp")
coarse_p

# save the plot
my_width <- 1
tex_width <- 5 * my_width; tex_height = 2.5 * my_width
save_tikz_plot(coarse_p, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-coarsening.tex")
