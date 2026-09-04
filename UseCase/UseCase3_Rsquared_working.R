# 1- Truth specification --------------------------------------------------
rm(list=ls())

pacman::p_load(future.apply, data.table, rethinking, dplyr, reticulate, mvtnorm)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)

#use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
#py_config() 
#source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
#source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")

{
  gen_AR1 <- function(p, rho = 0.6){
    # check if rank deficient: kappa value, SVD decomposition and eigenvalues not meeting tolerance level?
    
    Rxx <- matrix(NA,p,p)
    for(i in 1:p){
      for (j in 1:p){
        if(i == j) Rxx[i,j]<- 1
        if(i != j)  Rxx[i,j] <- rho^(abs(i-j))
      }
    }
    return(Rxx)
    
  }
  computeRsq <- function(M){
    # M is the correlation matrix of (Y,X)
    p = ncol(M)
    Rxy = as.matrix(M[2:p,1],ncol=1)
    Rxx = M[2:p, 2:p]
    
    Rsq = t(Rxy) %*% solve(Rxx) %*% Rxy
    return(Rsq)
  }
  
  make_mcmc<- function(S, asy.n, NCHAIN=1000, init= S, alpha=0.05, 
                       scale.init= 1, adaptive= TRUE, M= 100, acceptance.target = 0.44, 
                       seed= 123){
    p<- ncol(S)    #p: dimension of S (correlation matrix)
    d<- p*(p-1)/2  #d: dimension of correlation space
    accept <- 0
    q<- 0         # to trace scale_factor
    
    trace_scale<- array(NA_real_, dim=c(NCHAIN))
    scale_factor<- scale.init
    set.seed(seed)
    
    RCHAIN<- array(NA, dim= c(p,p,NCHAIN))
    RCHAIN[,,1] <- init
    
    for (i in 2:NCHAIN){
      
      Rcurrent<- RCHAIN[,,i-1]
      rcurrent<- Rcurrent[lower.tri(Rcurrent)]
      
      # Step 1: Draw R* (named as Rs)
      Psi_Rcurrent<- metaSEM::asyCov(Rcurrent, n = asy.n)
      
      Rs<- rwaldcloud(Rcurrent, Psi= scale_factor*Psi_Rcurrent, B = 1) %>% drop()
      rs<- Rs[lower.tri(Rs)]
      
      Psi_Rs<- metaSEM::asyCov(Rs, n= asy.n)
      # Step 2: Acceptance Ratio
      partA<- dmvnorm(x= rcurrent, mean = rs, sigma= Psi_Rs, log= T)
      partB<- dmvnorm(x= rs, mean= rcurrent, sigma= Psi_Rcurrent, log=T)
      logMH = partA - partB 
      
      partC<- wald.test(Rpop=Rs, Rsample= S, alpha= alpha, asy.n= asy.n)  # Rstar is at the center
      IR   <- partC$pval > alpha
      
      # Step 3: accept or reject
      if(IR==TRUE && logMH > log(runif(1))){ # accept Rs
        RCHAIN[,,i]<- Rs
        accept <- accept + 1
      } else {
        RCHAIN[,,i]<- Rcurrent
      }
      
      a_rate = accept/i
      cat("Acceptance Rate: ", round(a_rate,3),"\n")
      
      if(adaptive==TRUE){
        # update scale parameter every M steps
        ### If acceptance rate > 44% (acceptance.target parameter) 
        ############################ (step size is too small, then update scale to min(acceptance, 0.75)/ 0.44)
        ### If acceptance rate < 44% (step size is too big, then update scale to max(acceptance, 0.20)/ 0.44)
        ##### M: adaptive block
        if(i/M == floor(i/M)){
          a_rate=ifelse(a_rate>.75,.75,ifelse(a_rate<.1,.1,a_rate))
          scale_factor=scale_factor*a_rate/acceptance.target
          # cat("Scale Factor: ", scale_factor, "\n")
        }
      }
      trace_scale[i]<- scale_factor
    }
    attr(RCHAIN, "scale_trace") <- trace_scale
    
    return(RCHAIN)
  }
  
  simdriver<- function(Rxx, beta, k=0.85, ss=100, nchain = 1000, alpha = 0.05){
    
    p = ncol(Rxx)
    
    # Step 3: Compute S
    S = t(beta) %*% Rxx %*% (beta)
    
    # Step 4: Compute sigma2
    sigma2 = S *(1- k)/k
    
    # Step 5: Generate X ~ N(0, RXX) with sample size ss
    X <- rmvnorm(n = ss, mean = rep(0,p), sigma = Rxx)
    
    # Step 6: Generate Y = X'beta + eps where eps ~ N(0, sigma2)
    eps <- rnorm(n=ss, mean = 0, sd = sqrt(sigma2))
    Y <- X %*% (beta) + eps
    dt<- cbind(Y,X)
    
    # Step 7: estimate sample cor Rhat
    Rhat = cor(cbind(Y, X))
    
    # Step 8: map cor R to Rsq
    khat<- computeRsq(Rhat)
    
    # Step 9: Generate perturbation cloud
    RCHAIN <- make_mcmc(S= Rhat, asy.n= ss, NCHAIN = nchain, alpha = alpha)
    
    # Step 10: map RCHAIN to RSQ.hat.chain
    khat.chain <- sapply(1:nchain, function(i){
      computeRsq(RCHAIN[,,i])
    })
    
    khat.boot <- sapply(1:nchain, function(i){
      boot.id <-  sample(ss, ss, replace=TRUE)
      cor.boot<- cor(dt[boot.id,])
      computeRsq(cor.boot)
    })
    
    return(list(
      khat.point = khat,
      khat.chain = khat.chain,
      khat.boot = khat.boot
    ))
  }
  
  
}


# Step 0: Specify parameters
p<- 2
Rsq <- 0.6
ss <- 30
nchain = 10000


# Step 1: construct an AR1 Rxx corr structure 
Rxx<- gen_AR1(p=p, rho = 0.6)
Rxx
kappa(Rxx)  # should be fine as long as rho is not extremely large. Also set rho to large value so we can also study performance at the boundary 

# Step 2: Specify beta according to the signal scenario
beta <- c(-0.6, 0.6)
beta <- matrix(beta, nrow =p, ncol= 1)

# Step 3: Compute S
S = t(beta) %*% Rxx %*% (beta)

# Step 4: Compute sigma2
sigma2 = S *(1- Rsq)/Rsq

# Step 5: Generate X ~ N(0, RXX) with sample size ss
X <- rmvnorm(n = ss, mean = rep(0,p), sigma = Rxx)

# Step 6: Generate Y = X'beta + eps where eps ~ N(0, sigma2)
eps <- rnorm(n=ss, mean = 0, sd = sqrt(sigma2))
Y <- X %*% (beta) + eps

train.id<- sample(ss, round(ss*0.7,0))
test.id <- c(1:ss)[!(c(1:ss) %in% train.id)]



# Step 7: estimate sample cor Rhat
DT = cbind(Y,X)
Rhat = cor(DT)
Rxy = matrix(cor(Y, X), ncol = 1)
Rhat

Rhat.train <- cor(DT[train.id,])
Rhat.test <- cor(DT[test.id,])


##### Population R
{
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

# Step 8: map cor R to Rsq
computeRsq(Rhat)

# Step 9: Generate perturbation cloud

RCHAIN <- make_mcmc(S= Rhat, asy.n= ss, NCHAIN = nchain)
Rtrain.CHAIN <- make_mcmc(S= Rhat.train, asy.n= length(train.id), NCHAIN = nchain)
RCHAIN %>% dim

# Step 9B: Bootstrap based correlation (from step 6)
dt <- cbind(Y, X)

boot.id <-  sample(ss, ss, replace=TRUE)
boot.dt<- dt[boot.id,]
cor.boot<- cor(boot.dt)
computeRsq(cor.boot)

cor.boot <- sapply(1:nchain, function(i){
  boot.id <-  sample(ss, ss, replace=TRUE)
  cor.boot<- cor(dt[boot.id,])
  computeRsq(cor.boot)
})
summary( cor.boot) %>% round(.,3)
hist(cor.boot)
quantile(cor.boot, c(0.025, 0.975))



# Step 10: map RCHAIN to RSQ.hat.chain
Rsq.hat.chain <- sapply(1:nchain, function(i){
  computeRsq(RCHAIN[,,i])
})
Rsq.hat.chain %>% hist()
summary(Rsq.hat.chain)
range(Rsq.hat.chain) %>% diff


est.beta<- function(M){
  p <- ncol(M)-1
  Mxx<- M[2:(p+1), 2:(p+1)]
  Mxy<- M[1,2:(p+1)]
  beta<- solve(Mxx) %*% Mxy
  return(c(beta))
}
est.beta(Rhat.train)

computeRsq.test<- function(Mtest, beta.train){
  # beta.train is a row vector

  p <- ncol(Mtest)-1
  Mxx<- Mtest[2:(p+1), 2:(p+1)]
  Mxy<- Mtest[1,2:(p+1)]
  
  (2*beta.train) %*% Mxy - (beta.train %*% Mxx %*% matrix(beta.train, ncol=1))
}

Rsq.test.mcmc<- sapply(1:nchain, function(i){
  beta.i<- est.beta(Rtrain.CHAIN[,,i])
  computeRsq.test(Mtest = Rhat.test, beta.i)  
})

which.min(Rsq.test.mcmc)
summary(Rsq.test.mcmc)
hist(Rsq.test.mcmc)
Rtrain.CHAIN[,,i]

# Ensemble Simulation Driver
# Input: population correlation matrix, beta vector, true R squared, sample size ss
# Output: 1) point estimate of Rsq, 2) range estimates (eg 95% CI) of Rsq based on mcmc samples
# Notation: correlation matrix as M, Rsquared as k


p= 20
rho = 0.6
beta = sample(c(-0.4,0.4), p, replace= TRUE, prob = c(0.5,0.5))
ktrue = 0.8
ss= 100
nchain = 20000
alpha = 0.05

B = 300 # Repeat B times 
out.point <- array(NA, dim=c(1, B))
out.boot <- out.chain <- array(NA, dim=c(nchain, B))

M = gen_AR1(rho = rho, p= p)

out1<- simdriver(Rxx= M, beta= beta, k = ktrue, ss= ss, nchain = nchain, alpha= alpha)
out1$khat.point
out1$khat.chain %>% summary %>% round(.,3)
out1$khat.boot %>% quantile(., c(0.025, 0.975))

par(mfrow=c(2,1))
hist(out1$khat.chain, xlim= ktrue*c(0, 1.1), main="MCMC")
hist(out1$khat.boot, xlim= ktrue*c(0, 1.1), main="Bootstrap")


#for (i in 1:B){
#  out<- simdriver(Rxx= M, beta= beta, k = ktrue, ss= ss, nchain = nchain)
#  out.point[1,i] <- out$khat.point
#  out.chain[,i] <- out$khat.chain
#  out.boot[,i] <- out$khat.boot
#}

## Parallel run
plan(multisession, workers = parallelly::availableCores() - 1)
res_list <- future_lapply(1:B, function(i){
  library(mvtnorm)
  source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
  
  out<- simdriver(Rxx= M, beta= beta, k = ktrue, ss= ss, nchain = nchain)
  #out.point[1,i] <- out$khat.point
  #out.chain[,i] <- out$khat.chain
  #out.boot[,i] <- out$khat.boot
}, future.seed = TRUE)
plan(sequential)

# extract output
for (i in 1:B){
  out.point[1,i] <- res_list[[i]]$khat.point
  out.chain[,i] <- res_list[[i]]$khat.chain
  out.boot[,i] <- res_list[[i]]$khat.boot
}

# Report output
summary(c(out.point))

CI.mcmc<- apply(out.chain, 2, range)
CI.boot<- apply(out.boot, 2, function(x) quantile(x, c(0.025, 0.975)))
ktrue

mean((CI.mcmc[1,] < ktrue) & (CI.mcmc[2,] > ktrue))
mean((CI.boot[1,] < ktrue) & (CI.boot[2,] > ktrue))

out.chain %>% range
out.boot %>% quantile(.,c(0.025, 0.975))



# Scenario ----------------------------------------------------------------
### rho = 0.2, 0.5, 0.8. Use AR(1)
### for each rho setting, test all beta scenarios:
##### S1: all beta = 0.6
##### S2: all beta = 0.2
##### S3: mixed of 0.6 and 0.2
##### S4: one (or a few) beta = 0.8, rest = 0.2
##### S5: mixed of 0.4 and -0.4

nvar = 10
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

param_grid <- CJ(rho_id  = seq_along(rho_), beta_id = seq_along(beta_), ss_id   = seq_along(ss_) )
out.point <- array(NA_real_, dim = c(1, B, nrow(param_grid)))
out.boot  <- array(NA_real_, dim = c(nchain, B, nrow(param_grid)))
out.chain <- array(NA_real_, dim = c(nchain, B, nrow(param_grid)))


plan(multisession, workers = min(parallelly::availableCores() - 1, nrow(param_grid)))
for (j in seq_len(nrow(param_grid))) {
  rho_j  <- rho_[param_grid$rho_id[j]]
  beta_j <- beta_[[param_grid$beta_id[j]]]
  ss_j   <- ss_[param_grid$ss_id[j]]
  
  message("Running parameter set ", j, "/", nrow(param_grid))
  
  M_j <- gen_AR1(rho = rho_j, p = nvar)
  
  runs_j <- future_lapply(seq_len(B), function(b) {
    library(mvtnorm)
    source("docs/WorkingCode/Functions.R",local = TRUE)
    
    simdriver(Rxx = M_j, beta = beta_j, k = ktrue, ss = ss_j, nchain = nchain)
  }, future.seed = TRUE)
  
  for (i in seq_len(B)) {
    out.point[1, i, j] <- runs_j[[i]]$khat.point
    out.chain[, i, j] <- runs_j[[i]]$khat.chain
    out.boot[, i, j]  <- runs_j[[i]]$khat.boot
  }
  
  fout <- sprintf("sim_outputs/res_param_%03d.rds", j)
  
  saveRDS(
    list(
      j = j,
      param = param_grid[j, ],
      out.point = out.point[, , j, drop = FALSE],
      out.chain = out.chain[, , j, drop = FALSE],
      out.boot  = out.boot[, , j, drop = FALSE]
    ),
    file = fout
  )
  
  rm(runs_j)
  gc()
}
plan(sequential)

files <- list.files("sim_outputs", pattern = "\\.rds$", full.names = TRUE)
for (f in files) {
  x <- readRDS(f)
  j <- x$j
  out.point[, , j] <- x$out.point[, , 1]
  out.chain[, , j] <- x$out.chain[, , 1]
  out.boot[, , j]  <- x$out.boot[, , 1]
}


out.point[,,1]
out.chain[,,1]
out.boot[,,1]

CI.mcmc<- apply(out.chain, 2, range)
CI.boot<- apply(out.boot, 2, function(x) quantile(x, c(0.025, 0.975)))
ktrue

mean((CI.mcmc[1,] < ktrue) & (CI.mcmc[2,] > ktrue))
mean((CI.boot[1,] < ktrue) & (CI.boot[2,] > ktrue))

