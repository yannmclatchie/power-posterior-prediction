# set the seed
SEED <- 1234
set.seed(SEED)

# Define the true DGP
theta_ast <- 0.13

# repeat the experiment over different regimes
num_iters <- 1e3
iters <- 1:num_iters  # number of iterations
ns <- c(2, 10, 100)  # regimes of n
n_tau <- 1e2  # number of x-axis evaluations
taus <- c(0, 2^seq(-7, 7, length.out = n_tau), Inf)  # regimes of tau
priors <- list(#list(a = 1, b = 5, name = "weakly informative"),
               list(a = 1, b = 1, name = "uniform"))

# plotting
GR <- 1.61803
paper_theme <- (
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          legend.position="none"))
