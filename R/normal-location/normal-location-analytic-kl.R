library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
source("R/normal-location/config.R")

# compute the analytic KL risk
risk <- function(n, tau, sigma_0_2) {
  sigma_n_2 <- sigma_0_2 / (1 + n * tau * sigma_0_2)
  risk <- (1 / 2 * log(1 + sigma_n_2)
           + (1 + tau * sigma_n_2) / (2 * (1 + sigma_n_2))
           -1 / 2)
  return(list(risk = risk, n = n, tau = tau, sigma_0_2 = sigma_0_2))
}

# Evaluate the analytic risk across combinations
combis <- expand.grid(n = ns, tau = taus, sigma_0_2 = sigmas)
df <- combis |>
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
  geom_line(data = df,
            aes(tau, risk),
            colour = "black",
            size = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25) +
  facet_grid(sigma_0_2 ~ n, scales = "free_y",
             labeller = as_labeller(prior_names)) +
  scale_x_continuous(trans = "log10", breaks = 10^seq(-4, 4, length.out = 3)) +
  xlab("tau") +
  ylab("Risk") +
  paper_theme
p_risk

# save the plot
GR <- 1.61803
ggsave("./figs/normal-location-analytic-risk.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_risk, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-analytic-risk.tex")
