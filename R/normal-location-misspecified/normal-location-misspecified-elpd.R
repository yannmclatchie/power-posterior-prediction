library(simstudy)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
source("R/normal-location-misspecified/config.R")
source("R/normal-location-misspecified/utils.R")

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
n <- as.numeric(args[[1]]) # cycle through regimes of `n`
prior_name <- args[[2]] # cycle through prior
iter <- as.numeric(args[[3]]) # compute expectation with MC

# extract the prior from the dictionary
prior <- prior_dict[[prior_name]]

# iterate over values of tau
df <- taus|>
  map(\(tau) elpd_loo(iter = iter, 
                      n = n,
                      tau = tau,
                      prior = prior,
                      dgp = dgp)) |>
  bind_rows()

# save resutls to csv 
file_name <- paste0("data/normal-location-misspecified-elpd/",
                    "normal-location-misspecified-elpd-n_",
                    n,"-prior_",prior_name,"-iter_",iter,".csv")
write_csv(df, file = file_name)
