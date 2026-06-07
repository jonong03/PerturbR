library(data.table, dplyr)
files <- list.files("UseCase/msi/sim_outputs_v5", pattern = "\\.rds$", full.names = TRUE)
{
  # Param set
  library(data.table)
  nvar = 5
  ktrue = 0.9
  nchain = 20000
  B = 300
  alpha = 0.05
  
  rho_ = c(0.2, 0.5, 0.8)
  ss_ = c(200, 500, 1000)
  
  beta.s1 = rep(0.6, nvar)
  beta.s2 = rep(0.2, nvar)
  beta.s3 = sample(c(0.6, 0.2), nvar, replace=TRUE, prob = c(0.5,0.5))
  beta.s4 = beta.s2; beta.s4[sample(nvar,1)] <- 0.8
  beta.s5 = sample(c(-0.4,0.4), nvar, replace= TRUE, prob = c(0.5,0.5))
  beta_ = list(beta.s1, beta.s2, beta.s3, beta.s4, beta.s5)
  
  param_grid <- CJ( rho_id  = seq_along(rho_), beta_id = seq_along(beta_), ss_id   = seq_along(ss_) )
  
}


{
  out.point <- array(NA_real_, dim = c(1, B, length(files)))
  out.boot  <- array(NA_real_, dim = c(nchain, B, length(files)))
  out.chain <- array(NA_real_, dim = c(nchain, B, length(files)))
  for (f in files) {
    x <- readRDS(f)
    j <- x$j
  
    out.point[1, ,j] <- x$out.point  
    out.chain[, , j] <- x$out.chain
    out.boot[, , j]  <- x$out.boot
  }
  gc()
}

dim(out.point)


# Calibration summary output
apply(out.point,3, mean)
apply(out.point,3, sd)




param_grid[,Mean_PerturbationMean:= apply(apply(out.chain,c(2,3), mean),2,mean)]

out.chain[,,1]
dim(out.chain)
apply(apply(out.chain,c(2,3), mean),2,mean)

apply(out.chain,c(2,3), range) ## dim: 2, 300, 45
CI.mcmc <- apply(out.chain, c(2,3), range)
coverage.mcmc<- sapply(1:length(files), function(i) mean((CI.mcmc[1,,i] < ktrue) & (CI.mcmc[2,,i] > ktrue)))
coverage.mcmc
### what does the result imply?

CIwidth.mcmc<- sapply(1:length(files), function(i) mean(abs(CI.mcmc[1,,i]- CI.mcmc[2,,i])) )

CI.boot <- apply(out.boot, c(2,3), quantile, c(0.025, 0.975))
coverage.boot<- sapply(1:length(files), function(i) mean((CI.boot[1,,i] < ktrue) & (CI.boot[2,,i] > ktrue)))
coverage.boot
CIwidth.boot<- sapply(1:length(files), function(i) mean(abs(CI.boot[1,,i]- CI.boot[2,,i])) )

CIwidth.mcmc/ CIwidth.boot


param_grid[,TrueR2:= ktrue]
param_grid[,Mean_ObservedR2:= apply(out.point,3, mean)]
param_grid[,SD_ObservedR2:= apply(out.point,3, sd)]
param_grid[,bias:= Mean_ObservedR2 - TrueR2]
param_grid[,Mean_MonteCarloSD_ObservedR2:= apply(apply(out.chain, c(2,3),sd),2,mean)]
param_grid[,EmpCoverage_MCMC:= coverage.mcmc]
param_grid[,EmpCoverage_boot:= coverage.boot]
param_grid[,MeanWidth_MCMC:= CIwidth.mcmc]
param_grid[,MeanWidth_boot:= CIwidth.boot]
param_grid[,rel.width:= MeanWidth_MCMC/MeanWidth_boot]

round(param_grid, 3)

