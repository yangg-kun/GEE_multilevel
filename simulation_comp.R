## modeling by lm
theta_lm_all_seeds <- data.frame()
random_seed <- 0
for(i in 1:n_rep){
  random_seed <- random_seed + 1
  set.seed(random_seed)
  sdata <- simul_data(theta_r, 
                      n_smp, n_lv, 
                      corr = corr_coef)
  lm_outs <- lm_solve(sdata$V)
  theta_lm <- cbind(seed = random_seed, 
                    theta = "theta_lm", 
                    as.data.frame(t(unlist(lm_outs))))
  theta_lm_all_seeds <- rbind(theta_lm_all_seeds, theta_lm)
}


simul_lm <- theta_lm_all_seeds %>% 
  select(-c("seed", "theta")) %>% 
  colMeans() %>% 
  t() %>% as.data.frame() %>% `rownames<-`(c("theta_lm"))


## comparing GEE and lm in simulation
rename_gamma_names <- function(old_names){
  g0 <- which(startsWith(old_names, "gamma0"))
  g1 <- which(startsWith(old_names, "gamma1"))
  new_names <- old_names
  for(i in g0){
    new_names[i] <- 
      str_replace(new_names[i], "gamma0", "gamma0_")
  }
  for(i in g1){
    new_names[i] <- 
      str_replace(new_names[i], "gamma1", "gamma1_")
  }
  for(i in g1){
    new_names[i] <- 
      str_replace(new_names[i], "gamma1__inv", "gamma1_")
  }
  new_names
}
rename_gamma <- function(df){
  old_names <- colnames(df)
  new_names <- rename_gamma_names(old_names)
  colnames(df) <- new_names
  as.data.frame(df)
}


theta_all_seeds_ <- theta_all_seeds

del_names <- c("seed", "iter", "theta")
colnames(sd_all_seeds) <- 
  theta_all_seeds %>% 
  select(-del_names) %>% 
  colnames()

g1i <- startsWith(colnames(theta_all_seeds), "gamma1_inv")
theta_all_seeds_ <- theta_all_seeds %>% rename_gamma()
theta_all_seeds_[, g1i] <- 1 / theta_all_seeds_[, g1i]

theta_hat <- theta_all_seeds_ %>% 
  filter(theta == "theta_hat") %>% 
  select(-del_names) %>% 
  colMeans()
theta_real <- theta_all_seeds_ %>% 
  filter(theta == "theta_real") %>% 
  select(-del_names) %>% 
  colMeans()
simul_outs <- 
  rbind(theta_real, theta_hat)

sd_emp <- theta_all_seeds %>% 
  filter(theta == "theta_hat") %>% 
  select(-del_names) %>% 
  var() %>% diag() %>% sqrt()
sd_hat <- colMeans(sd_all_seeds) %>% 
  `names<-`(names(sd_emp))

simul_outs2 <- simul_lm %>% rename_gamma() %>% 
  bind_rows(as.data.frame(simul_outs), .) %>% 
  `rownames<-`(c("theta_real", 
                 "theta_gee", 
                 "theta_lm")) %>% 
  t()
simul_outs2[, -1] <- round(simul_outs2[, -1], digits = 3)

saving <- paste0("simul_outs.", corr_coef, ".csv")
simul_outs2 %>% 
  write.csv(saving, quote = F)


### boxplot
tmp <- theta_all_seeds_ %>% 
  filter(theta == "theta_hat") %>% 
  select(-c("seed", "iter", "theta")) %>% 
  select(-c("zeta", "eta")) %>% 
  gather(key = "theta", value = "estimate") %>% 
  mutate(theta = as.factor(theta)) %>% 
  mutate(method = "GEE")
tmp <- theta_lm_all_seeds %>% 
  filter(theta == "theta_lm") %>% 
  select(-c("seed", "theta")) %>% 
  select(-c("zeta", "eta")) %>% 
  gather(key = "theta", value = "estimate") %>% 
  mutate(theta = rename_gamma_names(theta)) %>% 
  mutate(theta = as.factor(theta)) %>% 
  mutate(method = "LM") %>% 
  rbind(tmp, .)
tmp_real <- rep(NA, nrow(tmp))
for(i in names(theta_real)){
  ri <- match(i, tmp$theta)
  if(is.na(ri)){
    next
  }
  # cat(i, theta_real[i], ri, "\n")
  tmp_real[ri] <- theta_real[i]
}
p <- tmp %>% 
  mutate(real = tmp_real) %>% 
  filter(str_starts(theta, "gamma")) %>% 
  mutate(theta = factor(theta, labels = c("gamma0", "gamma1"))) %>% 
  ggplot(aes(x = theta, y = estimate, 
             fill = method))
p + geom_boxplot(outlier.size = 0.8, 
                 outlier.alpha = 0.5) + 
  geom_point(aes(y = real), 
             color = "#B52129", 
             shape = 18, size = 3, 
             alpha = 0.8) + 
  my_theme + 
  ylim(c(-1.2, 2.4)) + 
  labs(x = paste0("CORR = ", corr_coef)) + 
  scale_fill_manual(values = c("#2171B5", "#E56D0B"))

saving <- paste0("simul_outs.", corr_coef, ".pdf")
ggsave(saving, width = 4, height = 4, dpi=128)