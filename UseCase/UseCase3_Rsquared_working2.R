# 1- Truth specification --------------------------------------------------
rm(list=ls())

pacman::p_load(xtable, future.apply, data.table, rethinking, dplyr, reticulate, mvtnorm)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
Rcpp::sourceCpp("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/calc_psi.cpp")

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
  datagen <- function(p, beta_value, ss, phi, R2true, ITER) {
    beta <- matrix(beta_value, nrow = p, ncol = 1)
    Rxx <- gen_AR1(p = p, rho = phi)
    S <- t(beta) %*% Rxx %*% beta
    sigma2 <- as.numeric(S) * (1 - R2true) / R2true
    
    n <- ss * ITER
    X <- rmvnorm(n=n, mean= rep(0,p), sigma= Rxx)
    eps <- rnorm(n=n, mean= 0, sd = sqrt(sigma2))
    
    Y <- X %*% beta + eps
    
    # Indices corresponding to each Monte Carlo replicate
    idx <- split(seq_len(n), rep(seq_len(ITER), each = ss))
    
    # Split into datasets
    DT.list <- lapply(
      idx,
      function(ii) as.matrix(cbind(Y = Y[ii], X[ii, , drop = FALSE]))
    )
    
    return(DT.list)
  }
  computeRsq_beta <- function(beta, theta) {
    2 * sum(beta * theta[-1, 1]) -
      drop(beta %*% theta[-1, -1] %*% beta)
  }
  split_DT <- function(DT.list, train_prop = 0.7) {
    
    lapply(DT.list, function(DT) {
      
      n <- nrow(DT)
      train_id <- sample.int(n, size = floor(train_prop * n))
      
      list(
        train = DT[train_id, , drop = FALSE],
        test  = DT[-train_id, , drop = FALSE]
      )
    })
  }
  estimate_beta <- function(theta) {
    theta_xx <- theta[-1, -1]
    theta_xy <- theta[-1, 1]
    beta_hat <- solve(theta_xx, theta_xy)
    as.numeric(beta_hat)
  }
}



# Parameter Matrix --------------------------------------------------------
{
  gc()
  R2true <- c(0.5)
  phi <- c(0.2, 0.5, 0.8)
  ss <- c(1000, 3000, 5000)
  p <- c(11, 14)
  
  beta.list <- c(
    #"function(p) rep(0.6, p)",
    #"function(p) rep(0.2, p)",
    "function(p) rep(c(-0.4, 0.4), length= p)"
  )
  
  param.m <- expand.grid(R2true = R2true, phi = phi, ss = ss, p = p, beta = beta.list,stringsAsFactors = FALSE) %>% as.data.table
  param.m$beta_value <- Map(
    function(f, p) {
      fun <- eval(parse(text = f))
      fun(p)
    },
    param.m$beta,
    param.m$p
  )
  param.m[, scenarioID := .I]
  #param.m <- param.m[rep(seq_len(nrow(param.m)), each = 3), ]
  #param.m[, method := rep(c("analytical", "bootstrap", "perturbR"), length.out = .N)]
}
{
  ncores <- detectCores()-1
  B = 500 # Number of Monte Carlo Replicates
  BOOT.ITER = 300   #J
  nchain = 3000
}
param.m


# Data Generation ---------------------------------------------------------
{
  param.m[, DT.list := Map(
    datagen,
    p, beta_value, ss, phi, R2true, ITER= B)]
}

# Analytical CI -----------------------------------------------------------

FL_eval <- function(DT.list, p, ss, R2true) {
  
  # Rhat = cor(DT), then compute R2
  R2hat <- sapply(DT.list, function(DT) computeRsq(cor(as.matrix(DT))))
  
  # CI for each Monte Carlo replicate
  fl.ci <- sapply(R2hat, function(r2) MBESS::ci.R2(R2 = r2, K = p, N = ss))
  lb <- unlist(fl.ci[1, ])
  ub <- unlist(fl.ci[3, ])
  
  c(coverage.fl = mean(R2true >= lb & R2true <= ub),
    width.fl    = mean(ub - lb))
}
fl.result <- Map(
  FL_eval,
  param.m$DT.list, param.m$p, param.m$ss, param.m$R2true)


# Bootstrap CI ------------------------------------------------------------
BOOT_eval <- function(DT.list, R2true, BOOT.ITER) {
  boot.ci <- lapply(DT.list, function(DT) {
    empR2dist <- sapply(seq_len(BOOT.ITER), function(i) {
      boot.id <- sample.int(nrow(DT), nrow(DT), replace = TRUE)
      DTboot <- DT[boot.id, , drop = FALSE]
      computeRsq(cor(DTboot))
    })
    quantile(empR2dist, c(0.025, 0.975))
  })
  
  boot.ci <- do.call(cbind, boot.ci)
  
  c(
    coverage.boot = mean(R2true >= boot.ci[1, ] & R2true <= boot.ci[2, ]),
    width.boot = mean(boot.ci[2, ] - boot.ci[1, ])
  )
}
boot.result <- parallel::mclapply(
  seq_len(nrow(param.m)),
  function(i) {
    tryCatch(
      BOOT_eval(DT.list = param.m$DT.list[[i]],R2true = param.m$R2true[i],BOOT.ITER = BOOT.ITER),
      error = function(e) {
        list(
          scenarioID = param.m$scenarioID[i],
          row = i,
          error = conditionMessage(e)
        )
      }
    )
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)
boot.result <- do.call(rbind, boot.result)
boot.result

# PerturbR ----------------------------------------------------------------
perturbR_eval <- function(DT.list, ss, R2true, nchain) {
  
  perturbR.ci <- lapply(DT.list, function(DT) {
    theta <- cor(DT)
    thetacloud <- make_mcmc_fast(S = theta, asy.n = ss, NCHAIN = nchain, alpha = 0.05)
    empR2dist <- vapply(seq_len(nchain), function(i) computeRsq(thetacloud[,,i]),numeric(1))
    quantile(empR2dist, c(0, 1))
  })
  
  perturbR.ci <- do.call(cbind, perturbR.ci)
  
  c(
    coverage.perturbR = mean(R2true >= perturbR.ci[1, ] & R2true <= perturbR.ci[2, ]),
    width.perturbR = mean(perturbR.ci[2, ] - perturbR.ci[1, ])
  )
}
perturbR.result <- parallel::mclapply(
  seq_len(nrow(param.m)),
  function(i) {
    
    tryCatch(
      perturbR_eval(
        DT.list = param.m$DT.list[[i]],
        ss = param.m$ss[i],
        R2true = param.m$R2true[i],
        nchain = nchain
      ),
      error = function(e) {
        list(
          scenarioID = param.m$scenarioID[i],
          row = i,
          error = conditionMessage(e)
        )
      }
    )
    
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)
perturbR.result <- do.call(rbind, perturbR.result)



# Assemble Results --------------------------------------------------------

out1<- copy(param.m[,.(scenarioID, R2true, phi, ss, p, beta_value)])
out1[,`:=`(
  coverage.fl = as.vector(sapply(fl.result, `[`, "coverage.fl")),
  width.fl    = as.vector(sapply(fl.result, `[`, "width.fl")),
  coverage.boot = boot.result[, "coverage.boot"],
  width.boot    = boot.result[, "width.boot"],
  coverage.perturbR = perturbR.result[, "coverage.perturbR"],
  width.perturbR    = perturbR.result[, "width.perturbR"]
)]
out1



#out1_p2p4<- copy(out1)
#out1_p9<- copy(out1)
out1_all<- rbind(out1_all, out1)
out1_highp <- rbind(out1, out1_highp)

saveRDS(out1_all, file="/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/UseCase/output/out1_all.rds")

setorder(out1_all, R2true, phi, p, ss)
out1_all[,beta_value:= sapply(beta_value, function(x) as.character(unique(x)))]
out1_all[,beta_value:= vapply(beta_value, paste, collapse=",", FUN.VALUE= character(1))]

xtable::xtable(out1_all)
out1_all[beta_value=="0.6" & phi ==0]

out1_all_xt<- copy(out1_all)
out1_all_xt<- out1_all_xt[,.(beta_value, phi, p, ss, coverage.fl, width.fl, coverage.boot, width.boot, coverage.perturbR, width.perturbR)]
out1_all_xt<- split(out1_all_xt, by="beta_value")
xtable.list <- lapply(out1_all_xt, function(dt) {xtable(dt[, !"beta_value"], digits = c(0, 2, 0, 0, 2, 2, 2, 2, 2, 2))})
xtable.list[[3]]


setorder(out1_highp, R2true, phi, p, ss)
out1_highp[phi==0]





library(ggplot2)

plotDT <- melt(
  out1_all,
  id.vars = c("scenarioID", "R2true", "phi", "ss", "p", "beta_value"),
  measure.vars = c("coverage.fl","coverage.boot","coverage.perturbR"),
  variable.name = "method",
  value.name = "coverage"
)

# plotDT<- plotDT[beta_value=="0.6" & phi==0.0]

plotDT[, method := factor(method,levels = c("coverage.fl", "coverage.boot", "coverage.perturbR"),labels = c("Analytical", "Bootstrap", "PerturbR"))]

# average across phi and beta
plot.avg <- plotDT[, .(coverage = mean(coverage)),by = .(ss, p, method)]

ggplot(
  plot.avg,
  aes(x = factor(ss), y = coverage, group = method, color = method)
  ) +
  geom_hline(yintercept = 0.95,linetype = "dashed")  +
  geom_line() +
  geom_point() +
  facet_wrap( ~ p , labeller = label_bquote(p == .(p))) +
  labs(x = "Sample size",y = "Empirical coverage",color = "Method") +
  theme_bw()















# Additional Analysis: OOS R2 ---------------------------------------------

# Split Data into train and test
param.m[, DT.split := lapply(DT.list, split_DT, train_prop = 0.7)]
param.m[, DT.train := lapply(DT.split, function(x) lapply(x, `[[`, "train"))]
param.m[, DT.test := lapply(DT.split, function(x) lapply(x, `[[`, "test"))]
param.m[, DT.split := NULL]

# Generate perturbations on DT.train and estimate theta_train
BOOT_eval_oos <- function(DT.train, DT.test, beta.true, BOOT.ITER) {
  out <- Map(function(train, test) {
    theta_test <- cor(test)
    # replicate-specific true OOS R2
    R2oos.true <- computeRsq_beta(beta = beta.true, theta = theta_test)
    
    # bootstrap training data only
    empR2oosdist <- vapply(seq_len(BOOT.ITER), function(b) {
      boot_id <- sample.int(nrow(train),nrow(train),replace = TRUE)
      theta_train <- cor(train[boot_id, , drop = FALSE])
      betahat_train <- estimate_beta(theta_train)
      computeRsq_beta(beta = betahat_train,theta = theta_test)
    }, numeric(1))
    
    boot.ci <- unname(quantile(empR2oosdist, c(0.025, 0.975)))
    
    c(
      coverage = as.numeric(R2oos.true >= boot.ci[1] &  R2oos.true <= boot.ci[2]),
      width = boot.ci[2] - boot.ci[1]
    )
    
  }, DT.train, DT.test)
  
  out <- do.call(rbind, out)
  
  c(
    coverage.boot.oos = mean(out[, "coverage"]),
    width.boot.oos    = mean(out[, "width"])
  )
}
boot.oos.result <- parallel::mclapply(
  seq_len(nrow(param.m)),
  function(i) {
    
    BOOT_eval_oos(
      DT.train = param.m$DT.train[[i]],
      DT.test = param.m$DT.test[[i]],
      beta.true = param.m$beta_value[[i]],
      BOOT.ITER = BOOT.ITER
      )
    
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)
boot.oos.result <- do.call(rbind, boot.oos.result)
boot.oos.result

perturbR_eval_oos <- function(DT.train, DT.test, beta.true, nchain) {
  
  out <- Map(function(train, test) {
    theta_train <- cor(train)
    theta_test  <- cor(test)
    
    # replicate-specific true OOS R2
    R2oos.true <- computeRsq_beta(beta = beta.true, theta = theta_test)
    
    # Perturb training correlation matrix
    thetacloud <- make_mcmc_fast(S = theta_train,asy.n = nrow(train),NCHAIN = nchain,alpha = 0.05)
    
    # theta_train uncertainty -> beta_hat -> OOS R2
    empR2oosdist <- vapply(seq_len(nchain), function(i) {
      beta_i <- estimate_beta(thetacloud[, , i])
      computeRsq_beta(beta = beta_i,theta = theta_test)
    }, numeric(1))
    perturbR.ci <- unname(quantile(empR2oosdist, c(0, 1)))
    
    c(
      coverage = as.numeric(R2oos.true >= perturbR.ci[1] &  R2oos.true <= perturbR.ci[2]),
      width = perturbR.ci[2] - perturbR.ci[1]
    )
  }, DT.train, DT.test)
  
  out <- do.call(rbind, out)
  c(
    coverage.perturbR.oos = mean(out[, "coverage"]),
    width.perturbR.oos    = mean(out[, "width"])
  )
}
perturbR.oos.result <- parallel::mclapply(
  seq_len(nrow(param.m)),
  function(i) {
    perturbR_eval_oos(DT.train = param.m$DT.train[[i]], DT.test = param.m$DT.test[[i]],beta.true = param.m$beta_value[[i]],nchain = nchain)
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)
perturbR.oos.result <- do.call(rbind, perturbR.oos.result)



# Assemble Results --------------------------------------------------------

out2<- copy(param.m[,.(scenarioID, R2true, phi, ss, p, beta_value)])
out2[, `:=`(
  coverage.boot.oos = boot.oos.result[, "coverage.boot.oos"],
  width.boot.oos    = boot.oos.result[, "width.boot.oos"]
)]
out2[, `:=`(
  coverage.perturbR.oos =
    perturbR.oos.result[, "coverage.perturbR.oos"],
  width.perturbR.oos =
    perturbR.oos.result[, "width.perturbR.oos"]
)]
out2















# Demo/ Visual Example
nchain= 10000
theta_hat = Rhat.list[[200]]

thetacloud<- make_mcmc_fast(S= theta_hat, asy.n = ss, NCHAIN= nchain, alpha = 0.05)
diff <- sapply(1:nchain, function(i){
  sqrt(sum((thetacloud[,,i]- Rhat.list[[1]])^2))
})
empR2dist <- sapply(1:nchain, function(i) computeRsq(thetacloud[,,i]))
bias = empR2dist - as.vector(computeRsq(theta_hat))
plot(diff, bias, ylab="R2_thetab - R2_thetahat", xlab = "Correlation Matrix Dispersion: ||theta_hat - theta_b||", main= "Diff in theta vs Diff in R2")
abline(h=0, col="red")


pacman::p_load(plotly)

R2hat1 = sapply(1:nchain, function(i) computeRsq(thetacloud[,,i]))
plot_ly(x=thetacloud[1,2,], y=thetacloud[1,3,], z=thetacloud[2,3,],
        type="scatter3d", mode="markers",
        marker=list(color=R2hat1, colorscale="Viridis", size=3, opacity=.4,
                   showscale=TRUE, colorbar=list(title="R2")), name="theta cloud") %>%
  add_markers(x=theta_hat[1,2], y=theta_hat[1,3], z=theta_hat[2,3],
              marker=list(color="red", size=7, symbol="diamond", showscale= FALSE), name="theta hat") %>%
  layout(scene=list(
    xaxis=list(title="theta12"),
    yaxis=list(title="theta13"),
    zaxis=list(title="theta23")
  ))








# Step 0: Specify parameters
p<- 3
Rsq <- 0.6
ss <- 30
nchain = 10000


# Step 1: construct an AR1 Rxx corr structure 
Rxx1<- gen_AR1(p=p, rho = 0.8)
Rxx2<- gen_AR1(p=p, rho = 0.1)
Rxx1
Rxx2
kappa(Rxx1)  # should be fine as long as rho is not extremely large. Also set rho to large value so we can also study performance at the boundary 
kappa(Rxx2)

# Step 2: Specify beta according to the signal scenario
beta <- c(-0.6, 0.6, -0.6)
beta <- matrix(beta, nrow =p, ncol= 1)

# Step 3: Compute S
S = t(beta) %*% Rxx2 %*% (beta)
S
# Step 4: Compute sigma2
sigma2 = S *(1- Rsq)/Rsq
sigma2

# Step 5: Generate X ~ N(0, RXX) with sample size ss
X <- rmvnorm(n = ss, mean = rep(0,p), sigma = Rxx2)

# Step 6: Generate Y = X'beta + eps where eps ~ N(0, sigma2)
eps <- rnorm(n=ss, mean = 0, sd = sqrt(sigma2))
Y <- X %*% (beta) + eps


train.id<- sample(ss, round(ss*0.7,0))
test.id <- c(1:ss)[!(c(1:ss) %in% train.id)]






pacman::p_load(MBESS)
MBESS::ci.R2(R2= computeRsq(Rhat), K=p, N=ss)
MBESS::ci.R2(R2= computeRsq(Rhat), K=p, N=ss)

MBESS::ci.R2(R2= Rsq, K=3, N=ss)


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

RCHAIN <- make_mcmc_fast(S= Rhat, asy.n= ss, NCHAIN = nchain)
Rtrain.CHAIN <- make_mcmc_fast(S= Rhat.train, asy.n= length(train.id), NCHAIN = nchain)
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




# Sim 1,2,3  ---------------------------------------------------------
R2true = 0.5
p=2
ss= 30
beta <- rep(c(0.6), length= p)
beta <- matrix(beta, nrow =p, ncol= 1)
sequential = FALSE # FALSE: Parallel run will be on
pacman::p_load(parallel)
ncores <- availableCores() - 1

ITER = 500
phi<- runif(ITER, 0.0, 0.0)
DT.list<- list()
Rhat.list <- list()
for (i in 1:ITER){
  Rxx<- gen_AR1(p=p, rho = phi[i])
  S = t(beta) %*% Rxx %*% (beta)
  sigma2 = S *(1- R2true)/R2true
  
  X <- rmvnorm(n = ss, mean = rep(0,p), sigma = Rxx)
  eps <- rnorm(n=ss, mean = 0, sd = sqrt(sigma2))
  Y <- X %*% (beta) + eps
  Y <- (Y - mean(Y))/ sd(Y)
  DT = cbind(Y,X)
  DT.list[[i]] <- DT
  
  Rhat = cor(DT)
  Rhat.list[[i]]<- Rhat
  
  # Compute population R2
  #Ryx<- (t(beta) %*% Rxx)/ as.vector(sqrt((t(beta) %*% Rxx %*% beta) + sigma2))
  #theta <- rbind(
  #  c(1, Ryx),
  #  cbind(t(Ryx), Rxx)
  #)
  #cat(computeRsq(theta),"\n")
}

## RES 1: FISHER/LEE

R2hat <- sapply(Rhat.list, computeRsq)
fl.ci <- sapply(R2hat, function(r2) MBESS::ci.R2(R2= r2, K=p, N=ss))
fl.ci.lb<- unlist(fl.ci[1,])
fl.ci.ub<- unlist(fl.ci[3,])

coverage.fl = mean(R2true >= fl.ci.lb & R2true <= fl.ci.ub )
width.fl = mean(fl.ci.ub  - fl.ci.lb)
coverage.fl
width.fl


## RES 2: BOOTSTRAP
BOOT.ITER= 3000
if(sequential == TRUE){
  boot.ci <- sapply(DT.list, function(DT){
    empR2dist<- sapply(1:BOOT.ITER, function(i){
      boot.id<- sample(ss, ss, replace= TRUE)
      DTboot <- DT[boot.id,]
      R2boot <- computeRsq(cor(DTboot))
      return(R2boot)
    })
    quantile(empR2dist, c(0.025, 0.975))
  })
  
}
boot.ci <- mclapply(DT.list, function(DT){
  empR2dist <- sapply(1:BOOT.ITER, function(i){
    boot.id <- sample(ss, ss, replace = TRUE)
    DTboot <- DT[boot.id, ]
    computeRsq(cor(DTboot))
  })
  
  quantile(empR2dist, c(0.025, 0.975))
}, mc.cores = ncores)
boot.ci <- do.call(cbind, boot.ci)

coverage.boot = mean(R2true >= boot.ci[1,] & R2true <= boot.ci[2,] )
width.boot = mean(boot.ci[2,] - boot.ci[1,])
coverage.boot
width.boot



## RES 3: PerturbR
nchain= 3000
if(sequential==TRUE){
  perturbR.ci <- sapply(Rhat.list, function(theta){
    thetacloud <- make_mcmc_fast(S= theta, asy.n = ss, NCHAIN = nchain, alpha = 0.05)
    empR2dist <- sapply(1:nchain, function(i) computeRsq(thetacloud[,,i]))
    quantile(empR2dist, c(0, 1))
  })
}

perturbR.ci <- mclapply(Rhat.list, function(theta){
  thetacloud <- make_mcmc_fast(
    S = theta, asy.n = ss,
    NCHAIN = nchain, alpha = 0.05
  )
  
  empR2dist <- sapply(1:nchain, function(i)
    computeRsq(thetacloud[,,i])
  )
  
  quantile(empR2dist, c(0, 1))
}, mc.cores = ncores)
perturbR.ci <- do.call(cbind, perturbR.ci)

coverage.perturbR = mean(R2true >= perturbR.ci[1,] & R2true <= perturbR.ci[2,] )
width.perturbR = mean(perturbR.ci[2,] - perturbR.ci[1,])
coverage.perturbR
width.perturbR

width.perturbR / width.fl

data.frame(
  Method=c("Analytical", "Bootstrap", "PerturbR"),
  Coverage= c(coverage.fl, coverage.boot, coverage.perturbR),
  Width = c(width.fl, width.boot, width.perturbR)
)




