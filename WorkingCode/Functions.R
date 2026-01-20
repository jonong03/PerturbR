# Functions
wald.test <- function(Rpop, Rsample, alpha = 0.05, asy.n = 100000, fisherz = FALSE) {
  # Rpop     : Target (true) correlation matrix
  # Rsample  : Estimated correlation matrix
  # alpha    : Type I error rate
  # indep    : If TRUE, use only marginal variances (diagonal Psi)
  # asy.n    : Sample size for asymptotic covariance
  # fisherz  : If TRUE, use Fisher-z transformation
  # pd_check : If TRUE, reject (return FALSE) when Rsample is not PD
  # pd_tol   : PD tolerance for smallest eigenvalue
  
  ## --- vectorize correlations ---
  Rpop_vec    <- Rpop[lower.tri(Rpop)]
  Rsample_vec <- Rsample[lower.tri(Rsample)]
  ps <- length(Rpop_vec)  # dimension = p*(p-1)/2
  
  ## --- Fisher z transform if requested ---
  if (fisherz) {
    center_vec <- atanh(Rpop_vec)
    sample_vec <- atanh(Rsample_vec)
    
    Psi <- asyCov.z(Rpop, asy.n = asy.n)     # asymptotic covariance in z-space
  } else {
    center_vec <- Rpop_vec
    sample_vec <- Rsample_vec
    
    Psi <- metaSEM::asyCov(Rpop, n = asy.n)  # asymptotic covariance in r-space
  }
  
  Psi0 <- Psi* asy.n   # covariance when sample size = 1, use this to study eigenstructure
  
  if(qr(Psi)$rank != ps) { return(c(distance=NA, df=NA, cricval=NA, sig=NA)) }
  
  ## --- Mahalanobis distance ---
  distance = mahalanobis(x= sample_vec,center = center_vec,cov = Psi)
  df= ps
  
  ## --- chi-square cutoff ---
  cricval<- qchisq(1 - alpha, df = df)
  pval<- as.numeric(pchisq(distance, df= df, lower.tail = F))
  reject= as.logical(distance >= cricval)
  out = c(T=distance, df= df, cricval= as.numeric(cricval), pval= pval, reject= reject) 
  
  ## --- ellipsoid membership ---
  return(out)
}
# Generate Sample (Full Spread)
rwaldcloud<- function(R0, n= 1e5, B=1, output="matrix"){
  r0 <- R0[lower.tri(R0)]
  Psi0<- metaSEM::asyCov(R0, n = 1)
  d<- ncol(R0)
  ps<- length(r0)
  
  # Eigen Decomposition and check if R0 (and Psi0) can be decomposed
  eig<- eigen(Psi0)
  if(!all(eig$values > 0)) {return("Covariance of R0 is not PSD")}
  U<- eig$vectors
  D<- diag(ps)
  diag(D)<- eig$values
  
  z<- mvtnorm::rmvnorm(n= B, sigma= diag(ps))  
  rhat<- r0 +U%*%sqrt(D/n)%*%matrix(z, ncol=B)
  
  Rhat<- array(NA, dim=c(d,d,B))
  for (b in 1:B){
    Rhatb <- diag(d)
    Rhatb[lower.tri(Rhatb)] <- rhat[,b]
    Rhatb[upper.tri(Rhatb)] <- t(Rhatb)[upper.tri(Rhatb)]
    Rhat[,,b]<- Rhatb
  }  
  
  # Check which is PSD (TRUE)
  is.psd<- sapply(1:B, function(b) all(eigen(Rhat[,,b])$values>0))
  non.psd.count<- sum(!is.psd)
  
  cat("non-PSD sample count:", non.psd.count, "\n")
  
  if(output=="matrix"){
    return(Rhat[,,is.psd])
  } else{
    return(t(rhat[,is.psd]))
  }
  
}
# Generate Sample (Uniform within boundary)
runifcloud<- function(R0, n= 1e5, B=1, alpha= 0.05, output="matrix"){
  r0 <- R0[lower.tri(R0)]
  Psi0<- metaSEM::asyCov(R0, n = 1)
  d<- ncol(R0)
  ps<- length(r0)
  
  # Eigen Decomposition and check if R0 (and Psi0) can be decomposed
  eig<- eigen(Psi0)
  if(!all(eig$values > 0)) {return("Covariance of R0 is not PSD")}
  U<- eig$vectors
  D<- diag(ps)
  diag(D)<- eig$values
  
  # rescaling factor
  q<- qchisq(p=1-alpha, df= ps)
  Z<- mvtnorm::rmvnorm(n= B, sigma= diag(ps)) 
  Zs<- Z/sqrt(rowSums(Z^2))   
  # check norm =1 : rowSums(Zs^2)=1
  r<- runif(B)^(1/ps) #random radius
  
  Rhat<- array(NA, dim=c(d,d,B))
  rhat_vech<- matrix(NA_real_, nrow = ps, ncol = B)
  
  for (b in 1:B){
    rhat<- r0 + U%*%sqrt(D/n)%*%matrix(Zs[b,],ncol=1)*(sqrt(q)*r[b])
    Rhatb <- diag(d)
    Rhatb[lower.tri(Rhatb)] <- rhat
    Rhatb[upper.tri(Rhatb)] <- t(Rhatb)[upper.tri(Rhatb)]
    
    Rhat[,,b]<- Rhatb
    rhat_vech[,b] <- rhat
  }  
  
  # Check which is PSD (TRUE)
  is.psd<- sapply(1:B, function(b) all(eigen(Rhat[,,b])$values>0))
  non.psd.count<- sum(!is.psd)
  
  cat("non-PSD sample count:", non.psd.count, "\n")
  
  if(output=="matrix"){
    return(Rhat[,,is.psd])
  } else{
    return(t(rhat_vech[,is.psd]))
  }
}