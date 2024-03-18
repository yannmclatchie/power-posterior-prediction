library(bayesflow)
source("R/normal-location-misspecified/config.R")

# central simulation script
simulate_data <- function(rep_id, n, mu_1, sigma_1, mu_2, sigma_2,
                          mu_3, sigma_3) {
  # define the DGP
  def <- defData(varname = "x1", formula = mu_1, variance = sigma_1)
  def <- defData(def, varname = "x2", formula = mu_2, variance = sigma_2)
  def <- defData(def, varname = "x3", formula = mu_3, variance = sigma_3)
  def <- defData(def, varname = "y", 
                 formula = "x1 | .3 + x2 | .4 + x3 | .3", 
                 dist = "mixture")
  
  # simulate from the DGP
  dx <- genData(n, def)
  y <- dx$y
  
  # return output
  return(list(rep_id = rep_id,
              n = n,
              y = y))
}

# initialise empty data list
datasets <- list()

# iterate over n
for (n in ns) {
  # print progress
  print(paste0("n = ",n))
  
  # simulate datasets
  dataset <- bayesflow::generate_from_dgp(dgp = simulate_data,
                                          n_datasets = num_iters, 
                                          n = n,
                                          mu_1 = mu_1,
                                          sigma_1 = sigma_1,
                                          mu_2 = mu_2,
                                          sigma_2 = sigma_2,
                                          mu_3 = mu_3,
                                          sigma_3 = sigma_3)
  
  # output to big list
  datasets[[paste(n)]] <- dataset
}

# save the datasets
saveRDS(datasets, "data/datasets/normal-mixture.RDS")  
