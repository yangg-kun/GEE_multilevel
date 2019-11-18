source("functions.R")

## parameters
n_smp <- 200  # sample size
n_lv <- 2  # number of levels(including standard level)
n_rep <- 1000  # times of repeats
corr_coef <- 0.2
random_seed <- 123
set.seed(random_seed)
theta_r <- rand_theta(n_lv)  # grenrating initial thetas

source("simulation_gee.R")
source("simulation_comp.R")