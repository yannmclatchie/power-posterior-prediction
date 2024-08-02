library(simstudy)
library(bayesflow)
library(tidyverse)
source("R/linear-regression/config.R")

# define the DGP 
simulate_data <- function(rep_id, n, p) {
  # generate linear regression data using simstudy
  def <- defRepeat(nVars = p, prefix = "x", formula = "0",
                   variance = "1", dist = "normal")
  def <- defData(def, varname = "yNormal", formula = lin_formula, 
                 dist = "normal", variance = "..sigma_ast")
  def <- defData(def, varname = "yOutlier", formula = 10, 
                 dist = "normal", variance = "..delta")
  def <- defData(def, varname = "y", dist = "mixture",
                 formula = "yOutlier | 0.2 + yNormal | 0.8") # NB: hard-coded eps
  dd <- genData(n, def)
  
  # convert to required format
  X <- as.matrix(dd)[, paste0("x", 1:p)]
  y <- dd$y
  
  # return output
  return(list(rep_id = rep_id,
              n = n, 
              p = p,
              X = X,
              y = y))
}

# initialise empty data list
datasets <- list()

# iterate over n
for (n in ns) {
  # print `n` regime
  print(paste0('n = ', n))
  
  # simulate datasets
  dataset <- bayesflow::generate_from_dgp(dgp = simulate_data,
                                          n_datasets = num_iters, 
                                          n = n,
                                          p = p)
  
  # output to big list
  datasets[[paste(n)]] <- dataset
}

# save the datasets
saveRDS(datasets, "data/datasets/linear_regression.RDS")  
