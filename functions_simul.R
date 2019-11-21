rand_d <- function(n = 1, min = 0, max = 1, digits = 0)
{
  round(runif(n, min, max), digits)
}

rand_theta <- function(n_lv)
{
  zeta <- rand_d(1, 0, 8, 1)
  eta <- rand_d(1, 0, 3, 1)
  gamma <- data.frame(gamma0 = 0, gamma1 = 1)
  gamma <- rbind(gamma, 
                 cbind(gamma0 = rand_d(n_lv - 1, -2, 2, 1), 
                       gamma1 = rand_d(n_lv - 1, 0, 2, 1))
  )
  theta <- list(zeta=zeta, eta=eta, 
                gamma0=gamma[[1]], 
                gamma1=gamma[[2]])
}

simul_data <- function(theta, n_smp=100, n_lv=2, corr=0.5)
{
  corr_V <- matrix(0, n_lv, n_lv)
  for(i in 1:n_lv){
    for(j in i:n_lv){
      if(i == j) next
      # corr_V[i, j] <- runif(1, 0, 1) # not positive definite
      corr_V[i, j] <- corr
    }
  }
  corr_V <- corr_V + t(corr_V) + diag(1, n_lv)
  
  var_V <- diag(theta$eta / (theta$gamma1 ^ 2))
  
  sigma_V <- sqrt(var_V) %*% corr_V %*% sqrt(var_V)
  mu_V <- (theta$zeta - theta$gamma0) / theta$gamma1
  
  V <- mvrnorm(n_smp, mu_V, sigma_V)
  list(V=V, theta=theta)
}
