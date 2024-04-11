library(bayesflow)
source("R/beta-binomial/config.R")

# central simulation script
simulate_data <- function(rep_id, n, theta_ast) {
  # simulate data from the true DGP
  y <- rbinom(n, 1, theta_ast)
  
  # return output
  return(list(rep_id = rep_id,
              n = n, 
              y = y,
              theta_ast = theta_ast))
}


# initialise empty data list
datasets <- list()

# iterate over n
for (n in ns) {
  # simulate datasets
  dataset <- bayesflow::generate_from_dgp(dgp = simulate_data,
                                          n_datasets = num_iters, 
                                          n = n,
                                          theta_ast = theta_ast)
  
  # output to big list
  datasets[[paste(n)]] <- dataset
}

# save the datasets
saveRDS(datasets, "data/datasets/beta_binomial.RDS")  
