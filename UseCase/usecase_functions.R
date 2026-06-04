
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
