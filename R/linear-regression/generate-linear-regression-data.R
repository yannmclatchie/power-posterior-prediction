library(simstudy)
library(bayesflow)
source("R/linear-regression/config.R")

# define the DGP 
simulate_data <- function(rep_id, n, p) {
  # generate linear regression data using simstudy
  def <- defRepeat(nVars = p, prefix = "x", formula = "0",
                   variance = "1", dist = "normal")
  def <- defData(def, varname = "y", formula = lin_formula, 
                 dist = "normal", variance = "..sigma_ast")
  def <- defData(def, varname = "isInlier", formula = "..eps",
                 dist = "binary")
  dd <- genData(n, def)
  
  # add in-liers in line with the Grunwald DGP
  dd <- dd |> 
    mutate_at("y", ~replace(., isInlier == 1, 0)) |>
    mutate(across(x1:x5, ~replace(., isInlier == 1, 0))) |>
    select(-isInlier)
  
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
