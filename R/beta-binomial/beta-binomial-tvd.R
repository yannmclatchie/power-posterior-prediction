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
  
  # compute the plug-in parameters
  theta_hat <- mean(y)

  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_beta_binom(theta_ast, alpha_n, beta_n)
  
  # compute the TVD triangle inequality components
  tvd_lemma <- tvd_beta_binom(theta_hat, alpha_n, beta_n)
  tvd_plug_in <- 0.5 * (abs(dbinom(0, 1, theta_ast)
                            - dbinom(0, 1, theta_hat))
                        + abs(dbinom(1, 1, theta_ast)
                              - dbinom(1, 1, theta_hat)))
  
  return(list(tvd = tvd,
              tvd_lemma = tvd_lemma,
              tvd_plug_in = tvd_plug_in,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# load in all datasets
datasets <- readRDS("data/datasets/beta_binomial.RDS")  

# restrict values of tau
taus <- taus[taus > 0 & is.finite(taus)]

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
  mutate(tvd = tvd * sqrt(n),
         tvd_lemma = tvd_lemma * sqrt(n),
         tvd_plug_in = tvd_plug_in * sqrt(n))

# Reduced number of iterations
df_100 <- df |>
  filter(iter <= 50)

# Produce ribbons for the figures
rdf <- df |>
  group_by(tau, n, prior) |>
  summarize(tvd_mean = mean(tvd),
            tvd_lemma_mean = mean(tvd_lemma),
            tvd_plug_in_mean = mean(tvd_plug_in),
            tvd_min = quantile(tvd, probs = 0.05),
            tvd_max = quantile(tvd, probs = 0.95))

# Plot the TVD over iterations
p_tvd <- ggplot() +
  geom_line(data = df_100,
            aes(tau, tvd, group = iter), 
            colour = "grey",
            alpha = 0.2) +
  geom_ribbon(data = rdf,
              aes(ymin = tvd_min,
                  ymax = tvd_max,
                  x = tau),
              colour = "black",
              alpha = 0.,
              linetype = "dotted") +
  facet_wrap( ~ n, scales = "fixed", nrow = 1) +
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

# make triangle inequality plotting df
triangle_df <- rdf |> 
  mutate(tvd_triangle = tvd_lemma_mean + tvd_plug_in_mean) |> # sanity check
  select(-c("tvd_min", "tvd_max")) |>
  reshape2::melt(id.vars = c("tau", "n", "prior"))

# replace underscores for plotting
triangle_df$variable <- gsub("_", " ", triangle_df$variable)

# plot the TVD components by triangle inequality
p_triangle <- triangle_df  |>
  ggplot(aes(tau, value, linetype = variable,
             alpha = variable,
             colour = variable)) +
  geom_hline(yintercept = 0, colour = "red") + 
  geom_line() + 
  facet_grid(prior ~ n, scales = "free_y") +
  scale_x_continuous(trans = "log2", 
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     label = function(x) ifelse(x == 0, "0", x)) + 
  scale_linetype_manual(values = c("longdash", "solid", "solid",
                                   "dotted")) +
  scale_alpha_manual(values = c(1, 1, 0.5, 1)) +
  scale_colour_manual(values = c("black", "black", "grey", "black")) +
  xlab("tau") +
  ylab("sqrt n TVD") +
  theme_classic() +
  theme(
    legend.background = element_rect(colour = "white", fill = "white"),
    legend.box.background = element_rect(colour = "white", fill = "white"),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 8))
p_triangle
