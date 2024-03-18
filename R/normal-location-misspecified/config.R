require("simstudy")
require("ggplot2")

# set the seed
SEED <- 1234
set.seed(SEED)

# define the true DGP
mu_1 <- -3
mu_2 <- 0
mu_3 <- 3
sigma_1 <- 1
sigma_2 <- 1
sigma_3 <- 1
lambda_1 <- 0.3
lambda_2 <- 0.4
lambda_3 <- 0.3

# compute the DGP sufficient statistics
mix_mean <- lambda_1 * mu_1 + lambda_2 * mu_2 + lambda_3 * mu_3
mix_mean2 <- (lambda_1 * (mu_1^2 + sigma_1) 
              + lambda_2 * (mu_2^2 + sigma_2) 
              + lambda_3 * (mu_3^2 + sigma_3))
mix_var <- mix_mean2 - mix_mean^2
sigma_mix <- sqrt(mix_var)

# repeat the experiment over different regimes
num_iters <- 1e4
iters <- 1:num_iters  # number of iterations
ns <- c(2, 10, 100)  # regimes of n
n_tau <- 1e2  # number of x-axis evaluations
taus <- 2^seq(-7, 7, length.out = n_tau)  # regimes of tau
# prior dictionary
weak_prior <- list(mu = 0, sigma = 1, name = "weak")
flat_prior <- list(mu = 0, sigma = sqrt(1e26), name = "flat")
prior_dict <- list("flat" = flat_prior, 
                   "weak" = weak_prior)

# plotting
GR <- 1.61803
paper_theme <- (
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          legend.position="none"))
