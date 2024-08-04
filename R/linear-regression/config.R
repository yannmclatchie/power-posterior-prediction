require("simstudy")
require("ggplot2")

# set the seed
SEED <- 1234
X_SEED <- SEED # for `x_sim`
set.seed(SEED)

# define the true DGP
p <- 5 # number of regression coefficients
# NB: below are _variances_
sigma_ast <- 1 # variance of bulk
delta <- 0.01 # variance of outlier
theta_ast <- c(rep(0.1, 4), rep(0, p - 4)) # true regression coefficients
lin_formula <- genFormula(theta_ast, sprintf("x%s", 1:p)) # bulk equation
contam_mu <- 0 # contamination location
eps <- 0.5 # contamination rate

# repeat the experiment over different regimes
num_iters <- 1e3
iters <- 1:num_iters  # number of iterations
ns <- c(50, 100, 500)  # regimes of n
n_tau <- 1e2  # number of x-axis evaluations
taus <- 2^seq(-7, 7, length.out = n_tau)  # regimes of tau
# prior dictionary
weak_prior <- list(mu_0 = 0, Sigma_0 = diag(p), name = "weak")
flat_prior <- list(mu_0 = 0, Sigma_0 = sqrt(1e26) * diag(p), name = "flat")
prior_dict <- list("flat" = flat_prior, 
                   "weak" = weak_prior)

# for Monte Carlo integration 
S <- 1e4

# plotting
GR <- 1.61803
paper_theme <- (
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          legend.position="none"))