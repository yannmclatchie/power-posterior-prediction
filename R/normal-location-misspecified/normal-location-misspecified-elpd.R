library(simstudy)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
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

# extract the data
datasets <- readRDS("data/datasets/normal-mixture.RDS")

# iterate over values of tau
df <- taus|>
  map(\(tau) elpd_loo(iter = iter, 
                      n = n,
                      tau = tau,
                      prior = prior,
                      datasets = datasets,
                      sigma_mix = sigma_mix)) |>
  bind_rows()

# save resutls to csv 
file_name <- paste0("data/normal-location-misspecified-elpd/",
                    "normal-location-misspecified-elpd-n_",
                    n,"-prior_",prior_name,"-iter_",iter,".csv")
write_csv(df, file = file_name)
