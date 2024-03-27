## ----
## KL divergence and LOO-CV elpd

elpd_loo <- function(iter, n, prior, tau, theta_ast, sigma_ast) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  X <- data$X
  y <- data$y
  
  # extract prior
  Sigma_0 <- prior$Sigma_0
  
  # compute the LOO-CV log predictive at each observation
  lpd <- numeric(n)
  for (i in 1:n) {
    lpd[i] <- lpd_loo_i(y, X, i, Sigma_0, sigma_ast, tau)
  }
  
  # compute the LOO-CV
  return(list(elpd = mean(lpd),
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}
lpd_loo_i <- function(y, X, i, Sigma_0, sigma_ast, tau) {
  # produce the deleted data
  y_loo <- y[-i]; X_loo <- as.matrix(X[-i,])
  y_oos <- y[i]; X_oos <- as.vector(X[i,])
  
  # fix dimensions in the n = 2 case
  if (length(y) == 2) {
    X_loo <- t(X_loo)
  }
  
  # compute the posterior parameters
  Sigma_n_inv <- tau / sigma_ast^2 * t(X_loo) %*% X_loo + solve(Sigma_0)
  Sigma_n <- solve(Sigma_n_inv)
  theta_n <- tau / sigma_ast^2 * Sigma_n %*% t(X_loo) %*% y_loo
  
  # compute the posterior predictive parameters
  mu_pred <- t(theta_n) %*% X_oos
  sigma_pred <- sqrt(t(X_oos) %*% Sigma_n %*% X_oos + sigma_ast^2)
  
  # evaluate the leave-one-out predictive
  lpd <- dnorm(x = y_oos, mean = mu_pred, sd = sigma_pred, log = TRUE)
  return(lpd)
}

## ----
## TVD

# simulation method for the predictors
x_sim <- function(S, p, eps) {
  # simulate multivariate normal
  x <- rmvnorm(S, mean = rep(0, p), sigma = diag(p))
  
  # return the predictors
  return(x)
}

# conditional log true DGP predictive density
log_mu_dens <- function(tilde_y, tilde_x, theta_ast, sigma_ast, eps, delta) {
  # compute the linear predictors
  y_pred <- tilde_x %*% theta_ast
  
  # compute and return the log predictive density
  lpd <- log(eps * dnorm(tilde_y, mean = y_pred, sd = sqrt(delta * sigma_ast^2)) 
             + (1 - eps) * dnorm(tilde_y, mean = y_pred, sd = sigma_ast))
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
                          sigma_ast, eps, delta) {
  mu_dens <- exp(log_mu_dens(tilde_y = y, tilde_x = x, theta_ast = theta_ast, 
                             sigma_ast = sigma_ast, eps = eps, delta = delta))
  pred_dens <- exp(log_pred_dens(tilde_y = y, tilde_x = x, theta_n = theta_n, 
                                 Sigma_n = Sigma_n, sigma_ast = sigma_ast))
  0.5 * abs(mu_dens - pred_dens)
}
# vectorise the integrand for dimensional coherence
tvd_integrand_vec <- Vectorize(tvd_integrand, vectorize.args = "y")
tvd_quad_int <- function(x, theta_n, Sigma_n, theta_ast, sigma_ast, eps, delta) {
  integrate(tvd_integrand_vec, x = x, theta_n = theta_n, Sigma_n = Sigma_n, 
            theta_ast = theta_ast, sigma_ast = sigma_ast, eps = eps, delta = delta,
            lower = -Inf, upper = Inf)$value
}

# define a loop-able function
tvd_quad_int_s <- function(s, tilde_x, theta_n, Sigma_n, theta_ast, 
                           sigma_ast, eps, delta){
  # extract simualted datum
  x <- tilde_x[s,]
  
  # evalute the TVD at this datum
  tvd_quad_int(x = x, 
               theta_n = theta_n,
               Sigma_n = Sigma_n, 
               theta_ast = theta_ast,
               sigma_ast = sigma_ast,
               eps = eps,
               delta = delta)
}

# main TVD experiment script
linreg_tvd <- function(iter, n, tau, prior, theta_ast, sigma_ast,
                       eps, delta) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  X <- data$X
  y <- data$y
  
  # extract prior
  Sigma_0 <- prior$Sigma_0
  
  # compute the posterior parameters
  Sigma_n_inv <- tau / sigma_ast^2 * t(X) %*% X + solve(Sigma_0)
  Sigma_n <- solve(Sigma_n_inv)
  theta_n <- tau / sigma_ast^2 * Sigma_n %*% t(X) %*% y
  
  # simulate the test predictors for MC integration
  tilde_x <- x_sim(S, p, eps)
  
  # compute the lpd at each test point
  MC_summands <- 1:S |> 
    purrr::map(\(s) tvd_quad_int_s(s = s,
                                   tilde_x = tilde_x, 
                                   theta_n = theta_n,
                                   Sigma_n = Sigma_n, 
                                   sigma_ast = sigma_ast,
                                   theta_ast = theta_ast,
                                   eps = eps,
                                   delta = delta))
  
  # evaluate the predictor-averaged TVD by Monte Carlo
  avg_tvd <- MC_summands |>
    unlist() |>
    mean() 
  
  # compute the TVD between the posterior predictive and the true DGP
  return(list(tvd = avg_tvd,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}
