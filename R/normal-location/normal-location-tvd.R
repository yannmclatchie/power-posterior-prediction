library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/normal-location/config.R")

# numerical computation of TVD between two normals
tvd_integrand <- function(x, mu_1, sigma_1, mu_2, sigma_2) {
  0.5 * abs(dnorm(x, mean = mu_1, sd = sigma_1)
            - dnorm(x, mean = mu_2, sd = sigma_2))
}
tvd_normal <- function(mu_1, sigma_1, mu_2, sigma_2) {
  integrate(tvd_integrand, mu_1, sigma_1, mu_2, sigma_2,
            lower = -Inf, upper = Inf)$value # rel.tol = 1e-9
}

# main experiment function
tau_tvd <- function(iter, n, prior, tau, theta_ast, sigma_ast) {
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
  
  # compute the plug-in parameters
  sigma_hat <- sigma_ast
  theta_hat <- bar_y
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_normal(theta_ast, sigma_ast, theta_post, sigma_pred)
  
  # compute the TVD triangle inequality components
  tvd_lemma <-tvd_normal(theta_hat, sigma_hat, theta_post, sigma_pred)
  tvd_plug_in <-tvd_normal(theta_ast, sigma_ast, theta_hat, sigma_ast)  
  
  return(list(tvd = tvd,
              tvd_lemma = tvd_lemma,
              tvd_plug_in = tvd_plug_in,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# load in all datasets
datasets <- readRDS("data/datasets/normal.RDS")  

# restrict values of tau
taus <- taus[taus > 0 & is.finite(taus)]

# evaluate the tvd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) tau_tvd(iter = iter, 
                                      n = n,
                                      tau = tau,
                                      prior = prior,
                                      theta_ast = theta_ast,
                                      sigma_ast = sigma_ast),
       .progress = TRUE) |>
  bind_rows()

# save resutls to csv 
file_name <- paste0("data/normal-location-tvd.csv")
#write_csv(df, file = file_name)

# read results
#df <- read_csv(file_name)

# fix ordering of priors
df$prior = factor(df$prior, levels=c('weak', 'flat'))

# scale the tvd by a function of n
df <- df |>
  mutate(tvd = tvd * sqrt(n),
         tvd_lemma = tvd_lemma * sqrt(n),
         tvd_plug_in = tvd_plug_in * sqrt(n))

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(tvd_mean = mean(tvd),
            tvd_lemma_mean = mean(tvd_lemma),
            tvd_plug_in_mean = mean(tvd_plug_in),
            tvd_min = quantile(tvd, probs = 0.05),
            tvd_max = quantile(tvd, probs = 0.95))

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Plot the TVD over iterations
p_tvd <- ggplot() +
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
  ylab("sqrt n TVD") +
  paper_theme
p_tvd

# save the plot
ggsave("./figs/normal-location-tvd.pdf", width = 5, height = 5 / GR)
my_width <- 1
tex_width <- 5 * my_width; tex_height = 2.5 * my_width
save_tikz_plot(p_tvd, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-tvd.tex")

# make triangle inequality plotting df
triangle_df <- rdf |> 
  mutate(tvd_triangle = tvd_lemma_mean + tvd_plug_in_mean) |> # sanity check
  select(-c("tvd_min", "tvd_max")) |>
  reshape2::melt(id.vars = c("tau", "n", "prior"))

# replace underscores for plotting
triangle_df$variable <- gsub("_", " ", triangle_df$variable)

# plot the TVD components by triangle inequality
p_triangle <- triangle_df  |>
  ggplot(aes(tau, value, linetype = variable,
             alpha = variable,
             colour = variable)) +
  geom_hline(yintercept = 0, colour = "red") + 
  geom_line() + 
  facet_grid(prior ~ n, scales = "free_y") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid",
                                   "dotted")) +
  scale_alpha_manual(values = c(1, 1, 0.5, 1)) +
  scale_colour_manual(values = c("black", "black", "grey", "black")) +
  xlab("tau") +
  ylab("sqrt n TVD") +
  theme_classic() +
  theme(#legend.position = c(0.86, 0.35),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.box.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 8))
p_triangle
