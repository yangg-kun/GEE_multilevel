---
GEE for multilevel data
---

This is a brief intorduction using GEE on simulation data as an example.

## generating simulation data

First, some parameters requires to be set.

```r
require(MASS)
source("functions_simul.R")  # loading required functions

n_smp <- 200  # sample size
n_lv <- 8  # number of levels (including base level)
# n_rep <- 1000  # times of repeats
corr_coef <- 0.2  # correlation coefficient for simulation data

## grenrating real values of thetas randomly
random_seed <- 123
set.seed(random_seed)
theta_r <- rand_theta(n_lv)
print(theta_r)
```

```
## $zeta
## [1] 2.3
## 
## $eta
## [1] 2.4
## 
## $gamma0
## [1]  0.0 -0.4  1.5  1.8 -1.8  0.1  1.6  0.2
## 
## $gamma1
## [1] 1.0 0.9 1.9 0.9 1.4 1.1 0.2 1.8
```

```r
sdata <- simul_data(theta_r, n_smp, n_lv, corr = corr_coef)
str(sdata)
```

```
## List of 2
##  $ V    : num [1:200, 1:8] 3.5527 1.367 -0.0718 -0.1507 5.2359 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : NULL
##  $ theta:List of 4
##   ..$ zeta  : num 2.3
##   ..$ eta   : num 2.4
##   ..$ gamma0: num [1:8] 0 -0.4 1.5 1.8 -1.8 0.1 1.6 0.2
##   ..$ gamma1: num [1:8] 1 0.9 1.9 0.9 1.4 1.1 0.2 1.8
```

`sdata` contains 2 parts:

* `V` is the data matrix for simulation($n_{sample}{\times}n_{level}$)  
* `theta` records ${\theta}s$  


note: `theta_r` is a list of 4 types of ${\theta}s$, so can also be determined manually as follow:

```r
## manually setting real values of thetas
theta_r <- list(zeta = 2.3, 
                gamma0 = c(0.0, -0.4,  1.5, 1.8, 
                           -1.8,  0.1, 1.6, 0.2), 
                eta = 2.4, 
                gamma1 = c(1.0, 0.9, 1.9, 0.9, 
                           1.4, 1.1, 0.2, 1.8))
```


## an example for modeling on simulation data with gee

`gee_solve` estimates ${\theta}s$ by Newton's method with moment estimating for initialization

* input:  
    + `V`: a matrix of raw data
    + `iter_times`: maximum times of iteration
    + `eps`: threshold for stopping iteration
* output:  
    + `theta_hat`: estimates of ${\theta}s$
    + `sd_hat`: estimates of standard deviates
    + `sigma_hat`: estimates of covariances
    + `theta_iters`: a list of estimates of ${\theta}s$ for each iteration


```r
M_raw <- sdata$V
n_lv <- ncol(M_raw)
source("functions_GEE.R")
esti_gee <- gee_solve(M_raw, iter_times = 1000, eps = 1e-6)
str(esti_gee, max.level = 2)
```

```
## List of 4
##  $ theta_iters:List of 3
##   ..$ :List of 4
##   ..$ :List of 4
##   ..$ :List of 4
##  $ theta_hat  :List of 4
##   ..$ zeta  : num 2.33
##   ..$ gamma0: num [1:7] -0.5221 1.5478 1.7089 -2.365 -0.0402 ...
##   ..$ eta   : num 2.73
##   ..$ gamma1: num [1:7] 0.956 2.027 1.012 1.56 1.187 ...
##  $ sigma_hat  : num [1:16, 1:16] 0.01363 0.01076 0.00978 0.01099 0.01112 ...
##   ..- attr(*, "dimnames")=List of 2
##  $ sd_hat     : num [1:16] 0.117 0.254 0.147 0.155 0.317 ...
```

`get_adj` transforms raw data to adjusted data with estimates of ${\theta}s$

```r
M_gee <- get_adj(M_raw, esti_gee$theta_hat)
```

`get_p` helps to get p-values

```r
p_raw <- get_p(m = M_raw, 
               test_f = friedman.test, 
               adjust = "BH", 
               n = n_lv - 1)
p_gee <- get_p(m = M_gee, 
               test_f = friedman.test, 
               adjust = "BH", 
               n = n_lv - 1)
sapply(list(raw = M_raw, adj = M_gee), 
       get_p, 
       test_f = friedman.test, 
       adjust = "BH", 
       n = n_lv - 1)
```

```
##           raw           adj 
## 1.403209e-101  1.000000e+00
```



## simulating for 1000 times

```r
n_smp <- 200  # sample size
n_lv <- 8  # number of levels (including base level)
n_rep <- 1000  # times of repeats
corr_coef <- 0.5  # correlation coefficient for simulation data

## grenrating real values of thetas randomly
random_seed <- 123
set.seed(random_seed)
theta_r <- rand_theta(n_lv)

theta_real <- 
  list(zeta = theta_r$zeta, 
       gamma0 = theta_r$gamma0[-1], 
       eta = theta_r$eta, 
       gamma1 = theta_r$gamma1[-1])
theta_real <- unlist(theta_real)

source("functions_GEE.R")

t0 <- proc.time()
theta_all_iterations <- c()
sd_all_iterations <- c()
for(i in 1:n_rep){
  set.seed(i)
  M_raw <- simul_data(theta_r, n_smp, n_lv, corr = corr_coef)
  M_raw <- M_raw$V
  esti_gee <- gee_solve(M_raw, eps = 1e-6)
  # M_gee <- get_adj(M_raw, esti_gee$theta_hat)
  theta_esti <- unlist(esti_gee$theta_hat)
  sd_esti <- esti_gee$sd_hat
  theta_all_iterations <- rbind(theta_all_iterations, theta_esti)
  sd_all_iterations <- rbind(sd_all_iterations, sd_esti)
}
t1 <- proc.time()
cat("time consuming:\n")
print(t1 - t0)
cat("\n")

theta_esti <- colMeans(theta_all_iterations)
sd_emp <- sqrt(sapply(as.data.frame(theta_all_iterations), var))
sd_esti <- colMeans(sd_all_iterations)
cbind(real = theta_real, 
      esti_gee = theta_esti, 
      sd_emp = sd_emp, 
      sd_gee = sd_esti)
```

```
## time consuming:
##    user  system elapsed 
##    2.06    0.01    2.07 
## 
##         real    esti_gee     sd_emp     sd_gee
## zeta     2.3  2.30068243 0.10571146 0.10930028
## gamma01 -0.4 -0.39680985 0.19731553 0.19575740
## gamma02  1.5  1.49961847 0.11605853 0.11934571
## gamma03  1.8  1.79885723 0.11151595 0.11304219
## gamma04 -1.8 -1.80582749 0.27312576 0.27179943
## gamma05  0.1  0.09176835 0.17570153 0.17212860
## gamma06  1.6  1.60112563 0.11473105 0.11709537
## gamma07  0.2  0.19594129 0.17265249 0.16761144
## eta      2.4  2.39535946 0.24106546 0.23759027
## gamma11  0.9  0.90025184 0.05605485 0.05458851
## gamma12  1.9  1.90449610 0.11891237 0.11600242
## gamma13  0.9  0.90242953 0.05699434 0.05469224
## gamma14  1.4  1.40215486 0.08809273 0.08516976
## gamma15  1.1  1.10430031 0.06670204 0.06688042
## gamma16  0.2  0.20056878 0.01226576 0.01216076
## gamma17  1.8  1.80453707 0.11232916 0.10952425
```

a comparison to linear regression

```r
lm_solve <- function(df){
  lm.outs <- list()
  for(i in 2:dim(df)[2])
  {
    colnames(df) <- NULL
    df <- as.data.frame(df)
    colnames(df)[1] <- "y"
    colnames(df)[i] <- "x"
    lm.out <- lm(y ~ x, df)
    lm.outs$zeta <- NA
    lm.outs$gamma0 <- c(lm.outs$gamma0, lm.out$coefficients[[1]])
    lm.outs$eta <- NA
    lm.outs$gamma1 <- c(lm.outs$gamma1, lm.out$coefficients[[2]])
  }
  lm.outs
}

t0 <- proc.time()
lm_all_iterations <- c()
for(i in 1:n_rep){
  set.seed(i)
  M_raw <- simul_data(theta_r, n_smp, n_lv, corr = corr_coef)
  M_raw <- M_raw$V
  esti_lm <- unlist(lm_solve(M_raw))
  lm_all_iterations <- rbind(lm_all_iterations, esti_lm)
}
t1 <- proc.time()
cat("time consuming:\n")
print(t1 - t0)
cat("\n")

theta_lm <- colMeans(lm_all_iterations)
cbind(real = theta_real, 
      esti_gee = theta_esti, 
      esti_lm = theta_lm)
```

```
## time consuming:
##    user  system elapsed 
##    5.48    0.00    5.48 
## 
##         real    esti_gee   esti_lm
## zeta     2.3  2.30068243        NA
## gamma01 -0.4 -0.39680985 0.9447903
## gamma02  1.5  1.49961847 1.8995063
## gamma03  1.8  1.79885723 2.0486025
## gamma04 -1.8 -1.80582749 0.2380667
## gamma05  0.1  0.09176835 1.1970727
## gamma06  1.6  1.60112563 1.9515307
## gamma07  0.2  0.19594129 1.2455036
## eta      2.4  2.39535946        NA
## gamma11  0.9  0.90025184 0.4524715
## gamma12  1.9  1.90449610 0.9539868
## gamma13  0.9  0.90242953 0.4530854
## gamma14  1.4  1.40215486 0.7042456
## gamma15  1.1  1.10430031 0.5517626
## gamma16  0.2  0.20056878 0.1000926
## gamma17  1.8  1.80453707 0.9048435
```

