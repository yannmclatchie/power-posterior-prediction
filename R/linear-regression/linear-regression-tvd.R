library(matrixStats)
library(mvtnorm)
source("R/linear-regression/config.R")

## ----
## Computation methods

# simulation method for the predictors
x_sim <- function(S, p, eps) {
  # simulate multivariate normal
  x <- rmvnorm(S, mean = rep(0, p), sigma = diag(p))
  
  # zero-inflation
  lambda <- runif(S)
  x[lambda < eps] <- 0
  
  # return the predictors
  return(x)
}

# conditional log true DGP predictive density
log_mu_dens <- function(tilde_y, tilde_x, theta_ast, sigma_ast, eps) {
  # compute the linear predictors
  y_pred <- tilde_x %*% theta_ast

  # compute and return the log predictive density
  if (tilde_y == 0) {
    lpd <- logSumExp(log(eps), 
                     log(1 - eps) + dnorm(0, mean = y_pred, 
                                          sd = sigma_ast, log = T))
  } else {
    lpd <- log(1 - eps) + dnorm(tilde_y, mean = y_pred, 
                                sd = sigma_ast, log = T)
  }
  return(lpd)
}

# conditional log posterior predictive density
log_pred_dens <- function(tilde_y, tilde_x, theta_n, Sigma_n, sigma_ast) {
  # compute the posterior predictive parameters
  y_pred <- tilde_x %*% theta_n
  sigma_pred <- sqrt(t(tilde_x) %*% Sigma_n %*% tilde_x + sigma_ast)

  # compute and return the log predictive density
  return(dnorm(tilde_y, mean = y_pred, sd = sigma_pred, log = T))
}

# TVD calculation by quadrature integration
tvd_integrand <- function(y, x, theta_n, Sigma_n, theta_ast, 
                          sigma_ast, eps) {
  0.5 * abs(log_mu_dens(tilde_y = y, tilde_x = x, theta_ast = theta_ast, 
                        sigma_ast = sigma_ast, eps = eps)
            - log_pred_dens(tilde_y = y, tilde_x = x, theta_n = theta_n, 
                            Sigma_n = Sigma_n, sigma_ast = sigma_ast))
}
# vectorise the integrand for dimensional coherence
tvd_integrand_vec <- Vectorize(tvd_integrand, vectorize.args = "y")
tvd_quad_int <- function(x, theta_n, Sigma_n, theta_ast, sigma_ast, eps) {
  integrate(tvd_integrand_vec, x = x, theta_n = theta_n, Sigma_n = Sigma_n, 
            theta_ast = theta_ast, sigma_ast = sigma_ast, eps = eps, 
            lower = -Inf, upper = Inf)$value
}

## ----
## Experiment

# use random parameters to test
theta_n <- rnorm(p)
toe <- toeplitz((p:1)/p)
Sigma_n <- rWishart(1, p, toe)[, , 1]

# simulate the test predictors
tilde_x <- x_sim(S, p, eps)

# try with one point
tvd_quad_int(x = tilde_x[1,], 
             theta_n = theta_n,
             Sigma_n = Sigma_n, 
             theta_ast = theta_ast,
             sigma_ast = sigma_ast,
             eps = eps)

# compute the lpd at each test point
tilde_x |> purrr::map(\(x) tvd_quad_int(x, 
                                        theta_n = theta_n,
                                        Sigma_n = Sigma_n, 
                                        sigma_ast = sigma_ast,
                                        theta_ast = theta_ast,
                                        eps = eps))

