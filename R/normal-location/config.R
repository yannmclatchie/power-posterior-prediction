# set the seed
SEED <- 1234
set.seed(SEED)

# define the true DGP
theta_ast <- 3
sigma_ast <- 1

# repeat the experiment over different regimes
num_iters <- 1e4
iters <- 1:num_iters  # number of iterations
ns <- c(2, 10, 100)  # regimes of n
n_lambda <- 1e2  # number of x-axis evaluations
taus <- 10^seq(-5, 5, length.out = n_lambda)
priors <- list(list(mu = 0, sigma = 1, name = "weak"),
               list(mu = 0, sigma = sqrt(1e26), name = "flat"))
sigmas <- c(1, 1e26)

# plotting
GR <- 1.61803
paper_theme <- (
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position="none"))
