library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bayesflow)
library(tidyverse)
source("R/linear-regression/config.R")
source("R/linear-regression/utils.R")

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
n <- as.numeric(args[[1]]) # cycle through regimes of `n`
prior_name <- args[[2]] # cycle through prior
iter <- as.numeric(args[[3]]) # compute expectation with MC

# extract the prior from the dictionary
prior <- prior_dict[[prior_name]]

# load in all datasets
datasets <- readRDS("data/datasets/linear_regression.RDS")

# Evaluate the elpd across combinations
df <- taus |>
  map(\(tau) elpd_loo(iter = iter, 
                      n = n,
                      tau = tau,
                      prior = prior,
                      theta_ast = theta_ast,
                      sigma_ast = sigma_ast),
      .progress = TRUE) |>
  bind_rows()
print("Done!")

# save resutls to csv 
file_name <- paste0("data/linear-regression-elpd/",
                    "linear-regression-elpd-n_",
                    n,"-prior_",prior_name,"-iter_",iter,".csv")
write_csv(df, file = file_name)
