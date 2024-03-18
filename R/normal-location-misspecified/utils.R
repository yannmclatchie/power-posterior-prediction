# numerical computation of TVD
dnormmix <- function(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                     sigma_3, lambda_1, lambda_2, lambda_3) {
  # compute the density of the gaussian mixture
  return(lambda_1 * dnorm(x, mean = mu_1, sd = sigma_1)
         + lambda_2 * dnorm(x, mean = mu_2, sd = sigma_2)
         + lambda_3 * dnorm(x, mean = mu_3, sd = sigma_3))
}
tvd_integrand <- function(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                          sigma_3, lambda_1, lambda_2, lambda_3,
                          mu_pred, sigma_pred) {
  0.5 * abs(dnorm(x, mean = mu_pred, sd = sigma_pred)
            - dnormmix(x, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                       sigma_3, lambda_1, lambda_2, lambda_3))
}
tvd_normal_misspecified <- function(mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                                    sigma_3, lambda_1, lambda_2, lambda_3,
                                    mu_pred, sigma_pred) {
  integrate(tvd_integrand, mu_1, sigma_1, mu_2, sigma_2, mu_3, 
            sigma_3, lambda_1, lambda_2, lambda_3,
            mu_pred, sigma_pred, lower = -Inf, upper = Inf)$value
}

# compute the TVD for a given tau/n/prior/dgp combo
tau_tvd <- function(iter, n, prior, tau, datasets, sigma_mix) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # compute the posterior parameters
  bar_y <- mean(y)
  sigma_post <- sqrt(1 / ((n * tau) / sigma_mix^2 + 1 / sigma_0^2))
  mu_pred <- sigma_post^2 * (mu_0 / sigma_0^2 + (n * tau * bar_y) / sigma_mix^2)
  sigma_pred <- sqrt(sigma_post^2 + sigma_mix^2)
  
  # compute the TVD between the posterior predictive and the true DGP
  tvd <- tvd_normal_misspecified(mu_1, sigma_1, mu_2, sigma_2, mu_3, 
                                 sigma_3, lambda_1, lambda_2, lambda_3,
                                 mu_pred, sigma_pred) 
  
  # return output
  return(list(tvd = tvd,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}

# compute the LOO-CV elpd for a given tau/n/prior/dgp combo
elpd_loo <- function(iter, n, prior, tau, datasets, sigma_mix) {
  # extract the data
  data <- datasets[[as.character(n)]][[iter]]
  y <- data$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # compute the LOO-CV log predictive at each observation
  lpd <- numeric(n)
  for (i in 1:n) {
    lpd[i] <- lpd_loo_i(y, i, mu_0, sigma_0, sigma_mix, tau)
  }
  
  # compute the LOO-CV elpd
  return(list(elpd = sum(lpd) / n,
              iter = iter,
              n = n,
              tau = tau,
              prior = prior$name))
}
lpd_loo_i <- function(y, i, mu_0, sigma_0, sigma_ast, tau) {
  # produce the deleted data
  y_loo <- y[-i]
  y_oos <- y[i]
  
  # compute posterior parameters
  n <- length(y)
  bar_y_i <- mean(y_loo)
  sigma_i <- sqrt(1 / (((n - 1) * tau) / sigma_ast^2 + 1 / sigma_0^2))
  mu_i <- sigma_i^2 * (mu_0 / sigma_0^2 + ((n - 1) * tau * bar_y_i) / sigma_ast^2)
  
  # evaluate the log predictive
  log_pred <- dnorm(x = y_oos, mean = mu_i, 
                    sd = sqrt(sigma_i^2 + sigma_ast^2), log = TRUE)
  return(log_pred)
}