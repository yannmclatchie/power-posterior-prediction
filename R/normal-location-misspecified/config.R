# set the seed
#SEED <- 1234
#set.seed(SEED)

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
# and produce a method for sampling from it
def <- defData(varname = "x1", formula = mu_1, variance = sigma_1)
def <- defData(def, varname = "x2", formula = mu_2, variance = sigma_2)
def <- defData(def, varname = "x3", formula = mu_3, variance = sigma_3)
def <- defData(def, varname = "y", 
               formula = "x1 | .3 + x2 | .4 + x3 | .3", 
               dist = "mixture")
mix_mean <- lambda_1 * mu_1 + lambda_2 * mu_2 + lambda_3 * mu_3
mix_mean2 <- (lambda_1 * (mu_1^2 + sigma_1) 
              + lambda_2 * (mu_2^2 + sigma_2) 
              + lambda_3 * (mu_3^2 + sigma_3))
mix_var <- mix_mean2 - mix_mean^2
dgp <- list(def = def, sigma_mix = sqrt(mix_var))

# repeat the experiment over different regimes
num_iters <- 1e4
iters <- 1:num_iters  # number of iterations
ns <- c(2, 10, 100)  # regimes of n
n_lambda <- 1e2  # number of x-axis evaluations
taus <- 10^seq(-5, 5, length.out = n_lambda)
#priors <- list(list(mu = 0, sigma = 1, name = "weak"),
#               list(mu = 0, sigma = sqrt(1e26), name = "flat"))

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
