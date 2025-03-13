library(bayesflow)
source("R/normal-location/config.R")

## Well-specified

# central simulation script
simulate_data <- function(rep_id, n, theta_ast, sigma_ast) {
  # simulate from DGP
  y <- rnorm(n, theta_ast, sigma_ast)
  
  # return output
  return(list(rep_id = rep_id,
              n = n, 
              y = y,
              theta_ast = theta_ast,
              sigma_ast = sigma_ast))
}

# initialise empty data list
datasets <- list()

# iterate over n
for (n in ns) {
  # simulate datasets
  dataset <- bayesflow::generate_from_dgp(dgp = simulate_data,
                                           n_datasets = num_iters, 
                                           n = n,
                                           theta_ast = theta_ast,
                                           sigma_ast = sigma_ast)
  
  # output to big list
  datasets[[paste(n)]] <- dataset
}

# save the datasets
saveRDS(datasets, "data/datasets/normal.RDS")  

## Misspecified

# central simulation script
simulate_bad_data <- function(rep_id, n, df_ast) {
  # simulate from DGP
  y <- rnorm(n, df_ast)
  
  # return output
  return(list(rep_id = rep_id,
              n = n, 
              y = y,
              df_ast = df_ast))
}

# initialise empty data list
datasets <- list()

# iterate over n
for (n in ns) {
  # simulate datasets
  dataset <- bayesflow::generate_from_dgp(dgp = simulate_bad_data,
                                          n_datasets = num_iters, 
                                          n = n,
                                          df_ast)
  
  # output to big list
  datasets[[paste(n)]] <- dataset
}

# save the datasets
saveRDS(datasets, "data/datasets/normal-misspecified.RDS")  
