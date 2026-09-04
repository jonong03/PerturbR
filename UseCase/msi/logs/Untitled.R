rm(list=ls()); gc()

library(future.apply)
library(data.table)
library(mvtnorm)
library(parallelly)
library(dplyr)

source("docs/WorkingCode/Functions.R")   # source required functions
source("UseCase/usecase_functions.R")      # source parameters for simulation
source("UseCase/msi/param_setup.R")      # source parameters for simulation



nvar < 5
rho_j <- 0.5
beta_j <- c(0.6, -0.6)
ss_j <- 200
nchain <- 15000
ktrue = 0.7
alpha= 0.05
B=50

M_j <- gen_AR1(rho = rho_j, p = nvar)

##### Population R
{
  Rxx <- M_j
  beta <- beta_j
  p<- nvar
  
  S = t(beta) %*% Rxx %*% (beta)
  sigma2 = S *(1- ktrue)/ktrue
  
  num = (t(beta) %*% Rxx)
  denom = sqrt( (t(beta) %*% Rxx %*% beta) + sigma2 )
  r_yx = as.vector(num)/ c(denom)
  Rpop = diag(p+1)
  Rpop[1,1] <- 1  #Ryy
  Rpop[2:(p+1), 2:(p+1)] <- Rxx
  Rpop[1,2:(p+1)] <- r_yx
  Rpop[2:(p+1), 1] <- r_yx
  Rpop
}



plan(sequential)
plan(multisession, workers = 6)
{
  
  runs_j <- lapply(1:B, function(j) simdriver(Rxx = M_j, beta = beta_j, k = ktrue, alpha= alpha, ss = ss_j, nchain = nchain, adaptive= F))
  
  Rsq.point <- array(NA_real_, dim = c(1, B))
  Rsq.train.mcmc<- Rsq.test.mcmc<- Rsq.mcmc <- array(NA_real_, dim = c(nchain, B))
  Rsq.boot  <- array(NA_real_, dim = c(nchain, B))
  
  for (i in seq_len(B)) {
    Rsq.point[1, i] <- runs_j[[i]]$khat.point
    Rsq.mcmc[, i] <- runs_j[[i]]$khat.chain
    Rsq.boot[, i]  <- runs_j[[i]]$khat.boot
    Rsq.train.mcmc[, i] <- runs_j[[i]]$Rsq.train.mcmc
    Rsq.test.mcmc[, i] <- runs_j[[i]]$Rsq.test.mcmc
    
  }
  
  Rhat <- runs_j[[1]]$Rhat
  Rmcmc <- runs_j[[1]]$mcmc
  Rboot <- runs_j[[1]]$boot

}
plan(sequential)
par(mfrow=c(1,2))
hist(Rsq.train.mcmc)
hist(Rsq.test.mcmc)
apply(Rsq.test.mcmc,2,summary) %>%round(.,3)

sapply(list(Rsq.mcmc, Rsq.boot, Rsq.train.mcmc, Rsq.test.mcmc), summary) %>% as.data.table

# Plotly ------------------------------------------------------------------
pacman::p_load(plotly)

# Layers:
# 1: centroid: True Theta
# 2: observation: Theta hat
# 3: Perturbations: Theta mcmc
library(rethinking)
fullspace = rlkjcorr(n=20000, K=3, eta=1)
plot_ly(x= Rpop[1,2], y= Rpop[1,3], z= Rpop[2,3], type = "scatter3d", mode="markers", color = I("black"), name="theta (") %>%
  add_markers(x= Rhat[1,2], y= Rhat[1,3], z= Rhat[2,3], type = "scatter3d", mode="markers", color = I("red"), name="theta hat") %>%
  add_markers(x= ~Rmcmc[1,2,], y= ~Rmcmc[1,3,], z= ~Rmcmc[2,3,], color= out.chain_j, colors="Viridis", name="theta mcmc", marker=list(opacity = 0.3, size=3)) %>%
  add_markers(x= ~Rboot[1,2,], y= ~Rboot[1,3,], z= ~Rboot[2,3,], color= out.boot_j, colors="Viridis", name="theta boot", marker=list(opacity = 0.3, size=3)) %>%
  add_markers(x= ~fullspace[,1,2], y= ~fullspace[,1,3], z= ~fullspace[,2,3], color= I("gray"), name="full space", marker=list(opacity = 0.2, size=2)) 

computeRsq(Rhat)
range(out.chain_j)
quantile(out.boot_j, c(0.025, 0.975))
Rpop

test.boot<- sapply(1:nchain, function(i) wald.test(Rsample= Rhat, Rpop = Rboot[,,i], asy.n = ss_j)) %>% t() %>% as.data.table
mean(test.boot$pval < 0.05)

test.mcmc<- sapply(1:nchain, function(i) wald.test(Rsample= Rhat, Rpop= Rmcmc[,,i], asy.n = ss_j)) %>% t() %>% as.data.table
mean(test.mcmc$pval < 0.05)


c(.701, 0.819) %*% solve(matrix(c(1,0.765, 0.765,1), nrow=2, ncol=2, byrow=T)) %*% matrix(c(0.701, 0.819), nrow=2, ncol=1)
