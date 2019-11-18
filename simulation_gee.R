## main process
theta_all_seeds <- data.frame()
sd_all_seeds <- data.frame()
random_seed <- 0

### uncomment to set default parameters
# theta_r <- list(zeta=6.4, gamma0=c(0, 1), 
#                 eta=0.2, gamma1=c(1, 0.4))

for(i in 1:n_rep){
  ### uncomment to monitor process
  # if((i %% (n_rep %/% 20)) == 0){
  #   cat(paste0(i, "/", n_rep, "\n"))
  # }
  random_seed <- random_seed + 1
  set.seed(random_seed)
  sdata <- simul_data(theta_r, n_smp, n_lv, corr = corr_coef)
  gee_seed <- 0
  while(T){
    gee_seed <- gee_seed + 1
    set.seed(gee_seed)
    gee_outs <- try(gee_solve(sdata$V, eps = 1e-3), silent = T)
    if(class(gee_outs) == "try-error"){
      next
    }
    break
  }
  # gee_outs <- gee_solve(sdata$V, eps = 1e-10)
  theta_real <- c(zeta = sdata$theta$zeta, 
                  gamma0 = sdata$theta$gamma0[-1], 
                  eta = sdata$theta$eta, 
                  gamma1_inv = 1 / sdata$theta$gamma1[-1])
  theta_real <- unlist(theta_real)
  theta_init <- unlist(gee_outs$theta[[1]])
  theta_hat <- unlist(gee_outs$theta[[length(gee_outs$theta)]])
  theta_all <- rbind(theta_real, theta_hat, theta_init)
  row.names(theta_all) <- NULL
  theta_all <- cbind(seed = random_seed, 
                     iter = length(gee_outs$theta) - 1, 
                     theta = c("theta_real", 
                               "theta_hat", 
                               "theta_init"), 
                     as.data.frame(theta_all))
  theta_all_seeds <- rbind(theta_all_seeds, theta_all)
  sd_all_seeds <- rbind(sd_all_seeds, gee_outs$sd)
  gee_outs$sigma
}

sd_emp <- theta_all_seeds %>% 
  filter(theta == "theta_hat") %>% 
  select(-c("seed", "iter", "theta")) %>% 
  var() %>% diag() %>% sqrt()
sd_hat <- colMeans(sd_all_seeds) %>% 
  `colnames<-`(colnames(sd_emp))
theta_hat <- theta_all_seeds %>% 
  filter(theta == "theta_hat") %>% 
  select(-c("seed", "iter", "theta")) %>% 
  colMeans()
theta_real <- theta_all_seeds %>% 
  filter(theta == "theta_real") %>% 
  select(-c("seed", "iter", "theta")) %>% 
  colMeans()

## saving estimates
simul_outs <- rbind(theta_real, theta_hat, sd_emp, sd_hat)
saving <- paste0("simul_gee.", corr_coef, ".csv")
write.csv(simul_outs, file = saving, quote = F)
