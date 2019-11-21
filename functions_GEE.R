init_theta <- function(v)
{
  var_v <- sapply(v, var)
  zeta_hat <- mean(v[, 1])
  eta_hat <- var(v[, 1])
  gamma1_hat <- sqrt(var_v[1] / var_v[-1])
  if(is.null(dim(v[, -1]))){
    gamma0_hat <- mean(v[, 1]) - mean(v[, -1]) * gamma1_hat
  }else{
    gamma0_hat <- mean(v[, 1]) - colMeans(v[, -1]) * gamma1_hat
  }
  theta_init <- list(zeta=zeta_hat, gamma0=gamma0_hat, 
                     eta=eta_hat, gamma1=gamma1_hat)
}

theta2list <- function(theta, theta_templ)
{
  theta_list <- list()
  j1 <- 0
  for(ni in names(theta_templ)){
    j0 <- j1 + 1
    j1 <- j0 - 1 + length(theta_templ[[ni]])
    theta_list[[ni]] <- theta[j0:j1]
  }
  theta_list
}

get_delta <- function(theta_, f, Vi)
{
  zeta_hat <- theta_$zeta
  gamma0_hat <- theta_$gamma0
  eta_hat <- theta_$eta
  gamma1_hat <- theta_$gamma1
  g <- c(c(1, 1 / gamma1_hat) * (zeta_hat - c(0, gamma0_hat)), 
         c(1, 1 / gamma1_hat) ^ 2 * (eta_hat + (zeta_hat - c(0, gamma0_hat)) ^ 2))
  Si <- f - matrix(1, ncol = 1, nrow = nrow(f)) %*% g
  Di_zeta <- c(1, 
               1 / gamma1_hat, 
               2 * zeta_hat, 
               2 * 1 / gamma1_hat ^ 2 * (zeta_hat - gamma0_hat))
  Di_gamma0 <- cbind(diag(c(1, 
                            -(1 / gamma1_hat))), 
                     diag(c(2 * zeta_hat, 
                            -2 * (zeta_hat - gamma0_hat) / (gamma1_hat ^ 2))))
  Di_eta <- c(0,
              0 * gamma1_hat, 
              1, 
              (1 / gamma1_hat) ^ 2)
  Di_gamma1 <- cbind(diag(c(0, 
                            -(zeta_hat - gamma0_hat) / (gamma1_hat ^ 2))), 
                     diag(c(0, 
                            -2 * (eta_hat + (zeta_hat - gamma0_hat) ^ 2) / (gamma1_hat ^ 3))))
  Di <- rbind(Di_zeta, Di_gamma0[-1,], Di_eta, Di_gamma1[-1,])
  U <- rowSums(Di %*% solve(Vi) %*% t(Si))
  W <- nrow(f) * Di %*% solve(Vi) %*% t(Di)
  delta <- solve(W) %*% U
  theta_hat <- unlist(theta_) + delta
  theta_hat <- theta2list(theta_hat, theta_)
  list(theta_init = theta_, f = f, Vi = Vi, 
       g = g, Si = Si, Di = Di, 
       U = U, W = W, delta = delta, 
       theta_hat = theta_hat)
}

get_sd <- function(delta_outs)
{
  Di <- delta_outs$Di
  Vi <- delta_outs$Vi
  Si <- delta_outs$Si
  n <- nrow(Si)
  B <- Di %*% solve(Vi) %*% t(Di)
  sigma <- t(solve(B)) %*% (
    (Di %*% solve(Vi) %*% t(Si) %*% Si %*% solve(Vi) %*% t(Di)) / n
  ) %*% solve(B) / n
  sd <- sqrt(diag(sigma))
  colnames(sigma) <- NULL
  rownames(sigma) <- NULL
  names(sd) <- NULL
  list(sigma = sigma, sd = sd)
}

gee_solve <- function(V, eps=1e-3, iter_times=30)
{
  theta_lists <- list()
  f <- as.matrix(cbind(V, V ^ 2))
  v <- as.data.frame(V)
  colnames(v) <- NULL
  var_v <- sapply(v, var)
  var2_v <- sapply(v, function(x){var(x ^ 2)})
  Vi <- diag(c(var_v, var2_v))
  theta_list <- init_theta(v)
  theta_lists[[1]] <- theta_list
  for(nt_i in 1:iter_times){
    delta_outs <- get_delta(theta_list, f, Vi)
    delta <- delta_outs$delta
    theta_list <- delta_outs$theta_hat
    theta_lists[[nt_i + 1]] <- theta_list
    if(t(delta) %*% delta < eps){
      break
    }
  }
  sd.theta <- get_sd(delta_outs)
  list(theta_iters = theta_lists, 
       theta_hat = theta_list, 
       sigma_hat = sd.theta$sigma, 
       sd_hat = sd.theta$sd)
}

get_adj <- function(df, theta_list){
  n_smp <- nrow(df)
  n_lv <- ncol(df)
  gamma1_hat <- rep(theta_list$gamma1, n_smp)
  gamma1_hat <- matrix(gamma1_hat, nrow = n_smp, byrow = T)
  gamma0_hat <- rep(theta_list$gamma0, n_smp)
  gamma0_hat <- matrix(gamma0_hat, nrow = n_smp, byrow = T)
  cbind(df[,1], df[,-1] * gamma1_hat + gamma0_hat)
}

get_p <- function(m, test_f, ..., adjust=NULL, n=NULL){
  colnames(m) <- NULL
  rownames(m) <- NULL
  m <- as.matrix(m)
  test_out <- test_f(m, ...)
  p <- test_out[["p.value"]]
  if(! is.null(adjust)){
    p <- p.adjust(p, adjust, n)
  }
  p
}
