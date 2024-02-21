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
tau_tvd <- function(iter, n, prior, tau, dgp) {
  # simulate data from the true DGP
  dx <- genData(n, dgp$def)
  y <- dx$y
  
  # extract prior
  mu_0 <- prior$mu
  sigma_0 <- prior$sigma
  
  # extract mixture variance
  sigma_mix <- dgp$sigma_mix
  
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