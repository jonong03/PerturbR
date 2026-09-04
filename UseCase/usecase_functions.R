
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
  est.beta<- function(M){
    p <- ncol(M)-1
    Mxx<- M[2:(p+1), 2:(p+1)]
    Mxy<- M[1,2:(p+1)]
    beta<- solve(Mxx) %*% Mxy
    return(c(beta))
  }
  computeRsq.test<- function(Mtest, beta.train){
    # beta.train is a row vector
    
    p <- ncol(Mtest)-1
    Mxx<- Mtest[2:(p+1), 2:(p+1)]
    Mxy<- Mtest[1,2:(p+1)]
    
    (2*beta.train) %*% Mxy - (beta.train %*% Mxx %*% matrix(beta.train, ncol=1))
  }
  
  simdriver<- function(Rxx, beta, k=0.85, ss=100, nchain = 1000, alpha = 0.05, adaptive = adaptive, acceptance.target = acceptance.target){
    
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
    RCHAIN <- make_mcmc(S= Rhat, asy.n= ss, NCHAIN = nchain, alpha = alpha, adaptive= adaptive, acceptance.target = acceptance.target)
    
    # Step 10: map RCHAIN to RSQ.hat.chain
    khat.chain <- sapply(1:nchain, function(i){
      computeRsq(RCHAIN[,,i])
    })
    
    RBOOT <- array(NA_real_, dim = c(p+1, p+1, nchain))
    for(i in 1:nchain){
      boot.id <-  sample(ss, ss, replace=TRUE)
      RBOOT[,,i]<- cor(dt[boot.id,])
    }
    khat.boot <- sapply(1:nchain, function(i){
      computeRsq(RBOOT[,,i])
    })
    
    ### For OOS R-square
    train.id<- sample(ss, round(ss*0.7,0))
    test.id <- c(1:ss)[!(c(1:ss) %in% train.id)]
    Rhat.train <- cor(dt[train.id,])
    Rhat.test <- cor(dt[-train.id,])
    Rhat.train.CHAIN <- make_mcmc(S= Rhat.train, asy.n= length(train.id), NCHAIN = nchain)
    Rsq.train.mcmc <- sapply(1:nchain, function(i){
      computeRsq(Rhat.train.CHAIN[,,i])
    })
    Rsq.test.mcmc<- sapply(1:nchain, function(i){
      beta.i<- est.beta(Rhat.train.CHAIN[,,i])
      computeRsq.test(Mtest = Rhat.test, beta.i)  
    })
    
    return(list(
      Rhat = Rhat,
      mcmc = RCHAIN,
      boot = RBOOT,
      khat.point = khat,
      khat.chain = khat.chain,
      khat.boot = khat.boot,
      Rsq.train.mcmc = Rsq.train.mcmc,
      Rsq.test.mcmc = Rsq.test.mcmc
    ))
  }
  
  
}
