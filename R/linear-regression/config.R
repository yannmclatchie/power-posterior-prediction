require("simstudy")
require("ggplot2")

# set the seed
SEED <- 1234
set.seed(SEED)

# define the true DGP
p <- 10
sigma_ast <- 1
theta_ast <- rep(0.1, p)
lin_formula <- genFormula(theta_ast, sprintf("x%s", 1:p))

# repeat the experiment over different regimes
num_iters <- 1e4
iters <- 1:num_iters  # number of iterations
ns <- c(2, 50, 200)  # regimes of n
n_tau <- 1e2  # number of x-axis evaluations
taus <- 2^seq(-7, 7, length.out = n_tau)  # regimes of tau
priors <- list(#list(mu_0 = 0, Sigma_0 = sqrt(1e26) * diag(p), name = "flat"),
               list(mu_0 = 0, Sigma_0 = diag(p), name = "weak"))

# plotting
GR <- 1.61803
paper_theme <- (
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          legend.position="none"))