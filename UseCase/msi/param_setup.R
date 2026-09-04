# Param set
library(data.table)
nvar = 2
ktrue = 0.5
nchain = 20000
B = 1
alpha = 0.05

rho_ = c(0.8)
ss_ = c(30, 50, 100)

beta.s1 = rep(0.6, nvar)
beta.s2 = rep(0.2, nvar)
beta.s3 = sample(c(0.6, 0.2), nvar, replace=TRUE, prob = c(0.5,0.5))
beta.s4 = beta.s2; beta.s4[sample(nvar,1)] <- 0.8
beta.s5 = sample(c(-0.4,0.4), nvar, replace= TRUE, prob = c(0.5,0.5))
beta_ = list(beta.s1, beta.s2, beta.s3, beta.s4, beta.s5)
beta_ = list(beta.s1)

param_grid <- CJ( rho_id  = seq_along(rho_), beta_id = seq_along(beta_), ss_id   = seq_along(ss_) )
