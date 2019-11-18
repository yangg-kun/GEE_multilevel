## functions basical
require(MASS)
require(tidyverse)


## functions for simulation
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


## functions for GEE
init_theta <- function(v)
{
  v <- as.matrix(v)
  var_v <- sapply(as.data.frame(v), var)
  zeta_hat <- mean(v[, 1])
  eta_hat <- var(v[, 1])
  gamma1_inv_hat <- sqrt(var_v[-1] / var_v[1])
  if(is.null(dim(v[, -1]))){
    gamma0_hat <- mean(v[, 1]) - mean(v[, -1]) / gamma1_inv_hat
  }else{
    gamma0_hat <- mean(v[, 1]) - colMeans(v[, -1]) / gamma1_inv_hat
  }
  theta_init <- list(zeta=zeta_hat, gamma0=gamma0_hat, 
                     eta=eta_hat, gamma1_inv=gamma1_inv_hat)
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
  gamma1_inv_hat <- theta_$gamma1_inv
  g <- c(c(1, gamma1_inv_hat) * (zeta_hat - c(0, gamma0_hat)), 
         c(1, gamma1_inv_hat) ^ 2 * (eta_hat + (zeta_hat - c(0, gamma0_hat)) ^ 2))
  Si <- f - matrix(1, ncol = 1, nrow = nrow(f)) %*% g
  Di_zeta <- c(1, 
               gamma1_inv_hat, 
               2 * zeta_hat, 
               2 * gamma1_inv_hat ^ 2 * (zeta_hat - gamma0_hat))
  Di_gamma0 <- cbind(diag(c(1, 
                            -gamma1_inv_hat)), 
                     diag(c(2 * zeta_hat, 
                            -2 * (zeta_hat - gamma0_hat) * gamma1_inv_hat ^ 2)))
  Di_eta <- c(0, 
              0 * gamma1_inv_hat, 
              1, 
              gamma1_inv_hat ^ 2)
  Di_gamma1 <- cbind(diag(c(0, 
                            zeta_hat - gamma0_hat)), 
                     diag(c(0, 
                            2 * (eta_hat + (zeta_hat - gamma0_hat) ^ 2) * gamma1_inv_hat)))
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
  names(sigma) <- NULL
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
  list(theta = theta_lists, 
       theta_hat = theta_list, 
       sigma = sd.theta$sigma, 
       sd = sd.theta$sd)
}

gee_test <- function(V, theta)
{
  n_smp <- nrow(V)
  n_lv <- ncol(V)
  gamma1 <- 1 / theta$gamma1_inv
  gamma1 <- rep(gamma1, n_smp)
  gamma1 <- matrix(gamma1, nrow = n_smp, byrow = T)
  gamma0 <- theta$gamma0
  gamma0 <- rep(gamma0, n_smp)
  gamma0 <- matrix(gamma0, nrow = n_smp, byrow = T)
  V_hat <- cbind(V[1], V[-1] * gamma1 + gamma0)
  test_out <- friedman.test(as.matrix(V_hat))
  list(V = V_hat, test_out=test_out)
}

get_adj <- function(df, theta_list){
  n_smp <- nrow(df)
  n_lv <- ncol(df)
  gamma1_hat <- theta_list$gamma1 %>% 
    rep(n_smp) %>% matrix(nrow = n_smp, byrow = T)
  gamma0_hat <- theta_list$gamma0 %>% 
    rep(n_smp) %>% matrix(nrow = n_smp, byrow = T)
  cbind(df[,1], df[,-1] * gamma1_hat + gamma0_hat)
}

get_p <- function(m, test_f, adjust=NULL, n=NULL){
  colnames(m) <- NULL
  rownames(m) <- NULL
  m <- as.matrix(m)
  test_out <- test_f(m)
  p <- test_out[["p.value"]]
  if(! is.null(adjust)){
    p <- p.adjust(p, adjust, n)
  }
  p
}


## functions for lm
lm_solve <- function(df){
  df <- sdata$V
  lm.outs <- list()
  lm.outs$zeta <- NA
  lm.outs$gamma0 <- c()
  lm.outs$eta <- NA
  lm.outs$gamma1 <- c()
  for(i in 2:dim(df)[2])
  {
    colnames(df) <- NULL
    df <- as.data.frame(df)
    colnames(df)[1] <- "y"
    colnames(df)[i] <- "x"
    lm.out <- lm(y ~ x, df)
    lm.outs$gamma0 <- c(lm.outs$gamma0, lm.out$coefficients[[1]])
    lm.outs$gamma1 <- c(lm.outs$gamma1, lm.out$coefficients[[2]])
  }
  lm.outs
}


## functions for plotting
require(ggsci)

my_theme <- 
  theme(text = element_text(size = 16, 
                            color = "black"), 
        plot.title = element_text(hjust = 0.2), 
        panel.grid =element_blank(), 
        panel.border = element_rect(fill = NA), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text = element_text(color = "black"), 
        # axis.text.x = element_text(angle = 15, vjust = 0.6), 
        # axis.title = element_blank(), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        legend.key = element_blank(), 
        # legend.position = "bottom", 
        # legend.direction = "horizontal", 
        # strip.text = element_blank(), 
        strip.background = element_rect(fill = "grey80", 
                                        color = "black")
  )

plot_den <- function(df){
  p <- df %>% as.data.frame() %>% 
    gather(key = "variable", value = "value") %>% 
    ggplot(aes(value, color = variable))
  p + stat_density(aes(color = variable, 
                       fill = variable), 
                   size = 1.2, alpha = 0.6, 
                   geom = "line", 
                   position = "identity") + 
    labs(x = NULL, y = NULL, color = "") + 
    guides(color = F) + 
    my_theme + 
    scale_color_lancet()
}
