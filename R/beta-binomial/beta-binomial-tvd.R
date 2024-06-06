library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(readr)
source("R/beta-binomial/config.R")

# numerical computation of TVD between the DGP and predictive
tvd_beta_binom <- function(p, alpha, beta) {
  # convert alpha and beta to m and s
  m <- alpha / (alpha + beta)
  s <- (alpha + beta)
  
  # compute the numerical summation
  0.5 * (abs(dbinom(0, 1, p)
             - dbetabinom(0, 1, m, s))
         + abs(dbinom(1, 1, p)
               - dbetabinom(1, 1, m, s)))
}

# main experiment function
tau_tvd <- function(iter, n, prior, tau, theta_ast) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y

  # extract prior
  alpha <- prior$a
  beta <- prior$b

  # compute posterior parameters
  x <- sum(y)
  z <- n - x
  alpha_n <- alpha + tau * x
  beta_n <- beta + tau * z
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_beta_binom(theta_ast, alpha_n, beta_n)
  return(list(tvd = tvd,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# load in all datasets
datasets <- readRDS("data/datasets/beta_binomial.RDS")  

# evaluate the elpd across combinations
combis <- expand.grid(iter = iters, n = ns, tau = taus, prior = priors)
df <- combis |>
  pmap(\(iter, n, tau, prior) tau_tvd(iter = iter, 
                                      n = n,
                                      tau = tau,
                                      prior = prior,
                                      theta_ast = theta_ast),
       .progress = TRUE) |>
  bind_rows()

# save results to csv 
file_name <- paste0("data/beta-binomial-tvd.csv")
#write_csv(df, file = file_name)

# read results from csv
df <- read_csv(file = file_name)

# scale the tvd by a function of n
df <- df |>
  mutate(tvd = tvd * sqrt(n))

# Group dataframe by iteration
gdf <- df |>
  group_by(n, tau, prior) |>
  summarize(tvd_mean = mean(tvd))

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(tvd_min = quantile(tvd, probs = 0.05),
            tvd_max = quantile(tvd, probs = 0.95))

# Produce vertical dashed lines
xdf <- df |>
  group_by(prior, n) |>
  summarise(xint = log(n) * sqrt((theta_ast * (1 - theta_ast)) / n)) |>
  distinct() 
ann_text <- data.frame(n = 2, 
                       prior = 'weak', 
                       lab = "Text",
                       tau = 1, tvd = 1)

# Plot the TVD over iterations
p_tvd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, tvd, group = iter), 
            colour = "grey",
            #size = 0.2, 
            alpha = 0.2) +
  geom_ribbon(data = rdf,
              aes(ymin = tvd_min,
                  ymax = tvd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              #size = 0.5,
              linetype = "dotted") +
  #geom_vline(data = xdf, aes(xintercept = xint), 
  #           linetype = "dashed") +
  #geom_text(data = ann_text, aes(x = tau, y = tvd), label = "sqrt(pq / n)") +
  facet_wrap( ~ n, scales = "fixed") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) +
  xlab("tau") +
  ylab("sqrt n TVD") +
  paper_theme
p_tvd

# save the plot
ggsave("./figs/beta-binomial-tvd.pdf", width = 5, height = 5 / GR)
my_width <- 0.8
tex_width <- 5 * my_width; tex_height = (5 / GR) * my_width
save_tikz_plot(p_tvd, width = tex_width, height = tex_height * 0.75,
               filename = "./tikz/beta-binomial-tvd.tex")
