# main experiment function
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
  
  # compute the TVD between the posterior predictive and the true DGP
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