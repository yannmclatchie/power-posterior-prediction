require("ggplot2")

# set the seed
SEED <- 12345
set.seed(SEED)

# define the true DGP
theta_ast <- 0
sigma_ast <- 1

# repeat the experiment over different regimes
num_iters <- 1e3
iters <- 1:num_iters  # number of iterations
ns <- c(2, 10, 100, 1000)  # regimes of n
n_tau <- 1e2  # number of x-axis evaluations
taus <- c(0, 2^seq(-7, 7, length.out = n_tau), Inf)  # regimes of tau
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
facet_names <- c(
  `elpd` = "elpd(tau)",
  `tau` = "log(tau)",
  `2` = "n = 2",
  `10` = "n = 10",
  `100` = "n = 100",
  `1000` = "n = 1000"
)
theme_sel <- theme_bw() +
  theme(axis.line.y=element_blank(),
        axis.line.x=element_line(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())