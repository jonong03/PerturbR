# Adaptive Metropolis-Hasting MCMC
# Target distribution: A set of correlation matrices that are true to observed S
#        In other words, which are the possible R0;s that could possibly generate S

# Parameters:
### NCHAIN
### S Matrix
### alpha: type1 error rate
### adaptive: indicator, TRUE or FALSE

# Output:
### RCHAIN
### acceptance rate (just print)
{
  rm(list=ls()); gc()
  pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr, plotly)
  
  use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
  py_config() 
  reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
  reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
  source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
  source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)
}

mcmc<- function(S, asy.n, NCHAIN=1000, init= S, alpha=0.05, 
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


{
  # Prep: Generate trueR
  p=3L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p3<- mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)

{
  # Prep: Generate trueR
  p=5L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p5<- mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)

{
  # Prep: Generate trueR
  p=5L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p10<- mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)

{
  # Prep: Generate trueR
  p=15L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p15_n500<- mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)
chain.p15_n50000<- mcmc(S=S, asy.n= asy.n*100, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)

{
  # Prep: Generate trueR
  p=20L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p20_n500<- mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)
chain.p20_n50000<- mcmc(S=S, asy.n= asy.n*100, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)

{
  # Prep: Generate trueR
  p=30L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p30_n500<- mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)
chain.p30_n50000<- mcmc(S=S, asy.n= asy.n*100, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)




# Compute reverse distance to S -------------------------------------------
compute_mdistance <- function(sample, S, asy.n, p = 0.95) {
  
  # Number of chains
  NCHAIN <- dim(sample)[3]
  
  # Precompute values from S (constant across chains)
  rs <- S[lower.tri(S)]
  Ps <- metaSEM::asyCov(S, n = asy.n)
  d  <- length(rs)
  c  <- qchisq(p = p, df = d)
  
  Qi <- sapply(1:NCHAIN, function(i) {
    
    rchain_mat <- sample[,,i]
    Pchain <- metaSEM::asyCov(rchain_mat, n = asy.n)
    rchain <- rchain_mat[lower.tri(rchain_mat)]
    
    # Wald distance
    D <- t(rchain - rs) %*% solve(Pchain) %*% (rchain - rs)
    
    # Final transformation
    u <- (D / c)^(d / 2)
    
    return(as.numeric(u))
  })
  
  return(Qi)
}

q.p20_n500<- compute_mdistance(chain.p20_n500, S=S, asy.n = 500)
q.p20_n50000<- compute_mdistance(chain.p20_n50000,  S=S, asy.n = 50000)
par(mfrow=c(2,2))
hist(q.p20_n500); hist(q.p20_n50000)
qqplot(x= q.p20_n500, y=runif(NCHAIN, 0,1)); qqline(y= q.p20_n500, distribution = function(p) qunif(p), col="red")
qqplot(x= q.p20_n50000, y=runif(NCHAIN, 0,1)); qqline(y= q.p20_n50000, distribution = function(p) qunif(p), col="red")

q.p15_n500<- compute_mdistance(chain.p15_n500)
q.p15_n50000<- compute_mdistance(chain.p15_n50000)




# Multiple Chain ----------------------------------------------------------
### Want to show that multiple chain with different starting points (or some parameters) will converge to the same stationary distribution
### Check trace plot: trace of each dimension should converge.

{
  # Prep: Generate trueR
  p=3L; ad=2L; N= 50000L; NCHAIN = 1e4
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  
  asy.n= 500L
  sampleid<- sample(N, asy.n)
  X <- Target$X[1,sampleid,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
chain.p3a<- mcmc(S=S, asy.n= asy.n, init=S, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)
chain.p3b<- mcmc(S=S, asy.n= asy.n, init=rlkjcorr(1,K= p, eta=1), NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)
chain.p3d<- mcmc(S=S, asy.n= asy.n, init=rlkjcorr(1,K= p, eta=1), NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)
chain.p3c<- mcmc(S=S, asy.n= asy.n, init=rlkjcorr(1,K= p, eta=1), NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)

library(parallel)

init_list <- list(
  p3a = S,
  p3b = rlkjcorr(1, K = p, eta = 1),
  p3d = rlkjcorr(1, K = p, eta = 1),
  p3c = rlkjcorr(1, K = p, eta = 1)
)

chain_list <- mclapply(
  init_list,
  function(init_val) {
    mcmc(
      S = S,
      asy.n = asy.n,
      init = init_val,
      NCHAIN = NCHAIN,
      acceptance.target = 0.4,
      M = 50
    )
  },
  mc.cores = 4
)

chain.p3a <- chain_list$p3a
chain.p3b <- chain_list$p3b
chain.p3d <- chain_list$p3d
chain.p3c <- chain_list$p3c


# Benchmark against bootstrap ---------------------------------------------
### Benchmark versus bootstrap on (1) Computing Cost & Time, (2) Performance Results
### model: 
pacman::p_load(microbenchmark, parallelly)

num_cores <- parallelly::availableCores(omit = 2)
my_cluster <- makeCluster(num_cores)
registerDoParallel(my_cluster)

res<- microbenchmark::microbenchmark(
  "Bootstrap n=50"= {boot.ind.seq(data1, ind_effect= the_ind_effect, B=B)},
  "Bootstrap n=500"= {boot.ind.seq(data2, ind_effect= the_ind_effect, B=B)},
  "Bootstrap n=5000"= {boot.ind.seq(data3, ind_effect= the_ind_effect, B=B)},
  "MCMC n=50"= {mcmc(S=S, asy.n= 50, init=S1, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)},
  "MCMC n=500"= {mcmc(S=S, asy.n= 500, init=S2, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)},
  "MCMC n=5000"= {mcmc(S=S, asy.n= 5000, init=S3, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 50)},
  times = 10) 
stopCluster(my_cluster)

res
autoplot(res)






q<- attr(chain1,"scale_trace")
plot(q, type="l")

### Autocorrelation plot
par(mfrow=c(2,3))
chain1<- chain.p20_n50000
acf(chain1[1,2,], lag.max = dim(chain1)[3])
acf(chain1[1,3,], lag.max = dim(chain1)[3])
acf(chain1[2,3,], lag.max = dim(chain1)[3])
acf(chain1[2,4,], lag.max = dim(chain1)[3])
acf(chain1[2,5,], lag.max = dim(chain1)[3])
acf(chain1[2,6,], lag.max = dim(chain1)[3])

### Trace plot
plot(chain1[1,2,], type="l")
plot(chain1[1,3,], type="l")
plot(chain1[1,4,], type="l")
plot(chain1[1,5,], type="l")
plot(chain1[1,6,], type="l")


X <- sapply(1:dim(chain1)[3], function(i) {
  chain1[,,i][lower.tri(chain1[,,i])]
})
X <- sapply(1:dim(chain2)[3], function(i) {
  chain2[,,i][lower.tri(chain2[,,i])]
})
library(plotly)
plot_ly(
  x = X[1,],
  y = X[2,],
  z = X[3,],
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 5,
    color = 1:dim(chain1)[3],              # <-- gradient by iteration
    colorscale = "Jet", # or "Plasma", "Jet"
    showscale = TRUE
  )
)

apply(X,1,median)[1:10] %>% round(.,2)
round(S[lower.tri(S)],2)[1:10]


# Check if samples are still uniform in adaptive setting
dim(chain1)
Qi <- sapply(1:NCHAIN, function(i){
  rchain<- chain1[,,i]
  Pchain <- metaSEM::asyCov(rchain, n= asy.n)
  rchain<- rchain[lower.tri(rchain)]
  
  rs <- S[lower.tri(S)]
  Ps <- metaSEM::asyCov(S, n= asy.n)
  
  D<- t(rchain - rs)%*%solve(Pchain)%*% (rchain-rs)
  d<- length(rs)
  c<- qchisq(p = 0.95, df = d)
  u<- (D/c)^(d/2)
  
  return(u)
})
hist(Qi)
qqplot(x= Qi, y=runif(NCHAIN, 0,1)); qqline(y= Qi, distribution = function(p) qunif(p), col="red")
?qqline
