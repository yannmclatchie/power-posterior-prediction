library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
source("R/normal-location/config.R")

# method to compute the Hellinger bounds
hel_bounds <- function(n, tau, prior, sigma_ast, theta_ast) {
  # compute the k0 factor
  k0 = sigma_ast^2 / prior$sigma^2
  
  # compute the bounds
  lb <- abs(sqrt((4 * sqrt(1 + 1 / (n * tau + k0))) 
                 / (2 * (2 + 1 / (n * tau + k0)) + n * tau^2 / (n * tau + k0)^2))
            * exp(-k0^2 * theta_ast^2 / 
                    (2 * sigma_ast^2 
                     * (2 
                        * (2 * (n * tau + k0) + 1) 
                        * (n * tau + k0) + n * tau^2)))
            - 2 / (sqrt(4 + 1 / n)))
  ub <- 2 * sigma_ast / sqrt(n * tau + k0)
  return(list(n = n, tau = tau, k0 = k0, ub = ub, lb = lb, prior = prior$name))
}

# regimes to compute the bounds over
combis <- expand.grid(n = ns, tau = taus, prior = priors)
hel_df <- combis |>
  pmap(\(n, tau, prior) hel_bounds(n = n,
                                   tau = tau,
                                   prior = prior,
                                   theta_ast = theta_ast,
                                   sigma_ast = sigma_ast)) |>
  bind_rows()

# consider only tau > 1
hel_df <- hel_df |> filter(tau > 1)

# Plot the Hellinger bounds
p_hel <- ggplot() +
  geom_ribbon(data = hel_df,
              aes(ymin = ub,
                  ymax = lb,
                  x = tau),
              colour = "grey",
              alpha = 0.1,
              linetype = "dotted") +
  facet_grid(prior ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log10", 
                     breaks = 10^seq(0, 4, length.out = 3)) +
  xlab("tau") +
  ylab("Hellinger") +
  paper_theme
p_hel

# save the plot
ggsave("./figs/normal-location-hel-bounds.pdf", width = 5, height = 5 / GR)
tex_width <- 5 * 0.8; tex_height = (5 / GR) * 0.8
save_tikz_plot(p_hel, width = tex_width, height = tex_height,
               filename = "./tikz/normal-location-hel-bounds.tex")
