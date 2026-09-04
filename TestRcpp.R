# Test Rcpp

rm(list=ls())

pacman::p_load(future.apply, data.table, rethinking, dplyr, reticulate, mvtnorm)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
Rcpp::sourceCpp("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/calc_psi.cpp")

exists("calc_psi_cpp")

R <- cor(rmvnorm(n=200, mean= rep(0, 10)))
Psi.R <- calc_psi(R, n = 300)
Psi.cpp <- calc_psi_cpp(R, n = 300)

max(abs(Psi.R - Psi.cpp))

microbenchmark::microbenchmark(
  R = calc_psi(R, 300),
  Rcpp = calc_psi_cpp(R, 300),
  times = 100
)


library(microbenchmark)
Rhat <- cor(rmvnorm(n=200, mean= rep(0, 2)))
ss = 300
nchain = 500
mb <- microbenchmark(
  original = make_mcmc(
    S = Rhat,
    asy.n = ss,
    NCHAIN = nchain,
    alpha = 0.05,
    seed = 123
  ),
  
  fast = make_mcmc_fast(
    S = Rhat,
    asy.n = ss,
    NCHAIN = nchain,
    alpha = 0.05,
    seed = 123
  ),
  
  times = 20
)

mb
summary(mb)
