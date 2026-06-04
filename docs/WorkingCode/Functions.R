# Functions

# Compute Covariance matrix for a correlation matrix R
calc_psi <- function(R, n = 1) {
  R <- as.matrix(R)
  p <- nrow(R)
  
  if (ncol(R) != p) { stop("R must be a square matrix.") }
  
  idx <- list()
  k <- 1
  
  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      idx[[k]] <- c(i, j)
      k <- k + 1
    }
  }
  
  k <- length(idx)
  psi <- diag(k)
  
  for (i in seq_len(k)) {
    e <- idx[[i]][1]
    f <- idx[[i]][2]
    
    psi[i, i] <- (1 - R[e, f]^2)^2
    
    if (i > 1) {
      for (j in seq_len(i - 1)) {
        g <- idx[[j]][1]
        h <- idx[[j]][2]
        
        tmp <- c(
          (R[e, g] - R[e, f] * R[f, g]) *
            (R[f, h] - R[f, g] * R[g, h]),
          
          (R[e, h] - R[e, g] * R[g, h]) *
            (R[f, g] - R[f, e] * R[e, g]),
          
          (R[e, g] - R[e, h] * R[h, g]) *
            (R[f, h] - R[f, e] * R[e, h]),
          
          (R[e, h] - R[e, f] * R[f, h]) *
            (R[f, g] - R[f, h] * R[h, g])
        )
        
        psi[i, j] <- 0.5 * sum(tmp)
        psi[j, i] <- psi[i, j]
      }
    }
  }
  
  return(psi/n)
}




wald.test <- function(Rpop, Rsample, alpha = 0.05, asy.n = 100000, fisherz = FALSE) {
  # Rpop     : Target (true) correlation matrix, this will be treated as the center and used to compute variance
  # Rsample  : Estimated correlation matrix
  # alpha    : Type I error rate
  # indep    : If TRUE, use only marginal variances (diagonal Psi)
  # asy.n    : Sample size for asymptotic covariance
  # fisherz  : If TRUE, use Fisher-z transformation
  # pd_check : If TRUE, reject (return FALSE) when Rsample is not PD
  # pd_tol   : PD tolerance for smallest eigenvalue
  
  asyCov.z <- function(R, asy.n= 10000) {
    # R: untransformed correlation matrix
    
    rvec <- R[lower.tri(R)]
    cov <- calc_psi(R, n= asy.n)
    
    #denominator: 1-rho^2
    #denominator[i,j] = (1-r_i^2)*(1-r_j^2)
    weights <- 1 - rvec^2
    denom_matrix <- outer(weights, weights)
    
    #Steiger defines psi = N * sigma
    psi_matrix <- asy.n * cov
    
    #Transformation: Eq10 and Eq11
    c_matrix <- psi_matrix / denom_matrix
    zcov_matrix <- c_matrix / (asy.n-3)
    diag(zcov_matrix)<- 1/ (asy.n-3)
    
    return(zcov_matrix)  
  }
  
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
    
    Psi <- calc_psi(Rpop, n = asy.n)  # asymptotic covariance in r-space
  }
  
  Psi0 <- Psi* asy.n   # covariance when sample size = 1, use this to study eigenstructure
  
  if(qr(Psi)$rank != ps) { return(c(T=NA, df=NA, cricval=NA, pval=NA, reject= NA)) }
  
  ## --- Mahalanobis distance ---
  distance = mahalanobis(x= sample_vec,center = center_vec,cov = Psi)
  df= ps
  
  ## --- chi-square cutoff ---
  cricval<- qchisq(1 - alpha, df = df)
  pval<- as.numeric(pchisq(distance, df= df, lower.tail = F))
  reject= as.logical(distance >= cricval)
  out = data.frame(T=distance, df= df, cricval= as.numeric(cricval), pval= pval, reject= reject) 
  
  ## --- ellipsoid membership ---
  return(out)
}
# Generate Sample (Full Spread)
rwaldcloud<- function(R0, n= 1e5, B=1, max.iter= B*10, output=c("full","lower"), scale=1){
  output <- match.arg(output)
  r0 <- R0[lower.tri(R0)]
  Psi0<- scale*calc_psi(R0, n = 1)
  d<- ncol(R0)
  ps<- length(r0)
  
  # Eigen Decomposition and check if R0 (and Psi0) can be decomposed
  eig<- eigen(Psi0)
  if(qr(Psi0)$rank != ps) {
    return(stop("Error: Covariance of R0 is not PSD"))
  }
  U<- eig$vectors
  D<- diag(ps)
  diag(D)<- eig$values
  
  #z<- mvtnorm::rmvnorm(n= B, sigma= diag(ps))  
  #rhat<- r0 +U%*%sqrt(D/n)%*%matrix(z, ncol=B)
  
  Rhat<- array(NA, dim=c(d,d,B))
  rhat.keep <- matrix(NA, nrow=ps, ncol=B)
  
  b <- 0
  iter <- 0
  non.psd.count <- 0
  
  while (b < B && iter < max.iter){
    iter <- iter + 1
    
    z<- mvtnorm::rmvnorm(n=1, sigma=diag(ps))
    rhatb<- r0 + U%*%sqrt(D/n)%*%matrix(z, ncol=1)
    
    Rhatb <- diag(d)
    Rhatb[lower.tri(Rhatb)] <- rhatb[,1]
    Rhatb[upper.tri(Rhatb)] <- t(Rhatb)[upper.tri(Rhatb)]
    
    if(all(eigen(Rhatb)$values > 0)){
      b <- b + 1
      Rhat[,,b] <- Rhatb
      rhat.keep[,b] <- rhatb[,1]
    } else{
      non.psd.count <- non.psd.count + 1
    }
  }
  
  cat("PSD sample count:", b, ";")
  cat("non-PSD sample count:", non.psd.count, "\n")
  
  if (b < B){
    stop(paste0("Could not generate ", B, " PSD samples within max.iter = ", max.iter, "."))
  }
  
  
  if(output=="full"){
    return(Rhat)
  } else{
    return(t(rhat.keep))
  }
  
}
rwaldcloud<- function(R0, Psi, n= 1e5, B=1, tol=1e-5, verbose=FALSE){
  # Psi must be defined at sample size = 1
  # If Psi is not defined, it will be estimated by calc_psi(R0, n= n)
  
  d<- ncol(R0)
  ps<- d*(d-1)/2
  out <- array(NA_real_, dim = c(d, d, B))
  
  if(missing(Psi)){
    Psi <- calc_psi(R0, n = n)
  }
  
  # Eigen Decomposition and check if R0 (and Psi0) can be decomposed
  eig<- eigen(Psi, symmetric = TRUE)
  if (any(eig$values < -tol)) {
    stop("Psi0 is not positive semidefinite.")
  }
  
  A <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = ps)
  
  r0 <- R0[lower.tri(R0)]
  b <- 0
  iter <- 0
  non.psd.count <- 0
  max.iter <- 100
  
  while (b < B && iter < max.iter) {
    iter <- iter + 1
    
    z <- matrix(rnorm(ps), ncol = 1)
    rhatb <- as.vector(r0 + A %*% z)
    
    Rhatb <- diag(d)
    Rhatb[lower.tri(Rhatb)] <- rhatb
    Rhatb[upper.tri(Rhatb)] <- t(Rhatb)[upper.tri(Rhatb)]
    
    evals <- eigen(Rhatb, symmetric = TRUE, only.values = TRUE)$values
    
    if (all(evals > tol)) {
      b <- b + 1
      out[, , b] <- Rhatb
    } else {
      non.psd.count <- non.psd.count + 1
    }
  }
  
  if (verbose) {
    cat("Accepted:", b, "; Rejected:", non.psd.count,
        "; Total iterations:", iter, "\n")
  }
  
  if (b < B) {
    stop(sprintf("Could not generate %d positive-definite samples within max.iter = %d.", B, max.iter))
  }
  
  #attr(out, "accepted") <- b
  #attr(out, "rejected") <- non.psd.count
  #attr(out, "iterations") <- iter
  
  return(out)
}

# Generate Sample Randomly (1-alpha spread)
runifcloud<- function(R0, n= 1e5, B=1, alpha= 0.05, output=c("full","lower")){
  output <- match.arg(output)
  r0 <- R0[lower.tri(R0)]
  Psi0<- calc_psi(R0, n = 1)
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
  radius <- runif(B)
  r<- radius^(1/ps) #random radius
  
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
# Generate Sample Randomly (1-alpha spread) # updated
runifcloud<- function(R0, n= 1e5, B=1, alpha= 0.05, max.iter= B*10, output=c("full","lower")){
  output <- match.arg(output)
  r0 <- R0[lower.tri(R0)]
  Psi0<- calc_psi(R0, n = 1)
  d<- ncol(R0)
  ps<- length(r0)
  
  # Eigen Decomposition and check if R0 (and Psi0) can be decomposed
  eig<- eigen(Psi0)
  if(qr(Psi0)$rank != ps) {
    return(stop("Error: Covariance of R0 is not PSD"))
  }
  U<- eig$vectors
  D<- diag(ps)
  diag(D)<- eig$values
  
  # rescaling factor
  q<- qchisq(p=1-alpha, df= ps)
  
  Rhat<- array(NA, dim=c(d,d,B))
  rhat.keep <- matrix(NA, nrow=ps, ncol=B)
  
  b <- 0
  iter <- 0
  non.psd.count <- 0
  
  while (b < B && iter < max.iter){
    iter <- iter + 1
    
    Z<- mvtnorm::rmvnorm(n=1, sigma=diag(ps))
    Zs<- as.numeric(Z / sqrt(sum(Z^2)))   # unit direction # check norm =1 : rowSums(Zs^2)=1
    radius <- runif(1)
    r <- radius^(1/ps)                    # random radius
    
    rhatb<- r0 + U%*%sqrt(D/n)%*%matrix(Zs, ncol=1)*(sqrt(q)*r)
    
    Rhatb <- diag(d)
    Rhatb[lower.tri(Rhatb)] <- rhatb[,1]
    Rhatb[upper.tri(Rhatb)] <- t(Rhatb)[upper.tri(Rhatb)]
    
    if(all(eigen(Rhatb)$values > 0)){
      b <- b + 1
      Rhat[,,b] <- Rhatb
      rhat.keep[,b] <- rhatb[,1]
    } else{
      non.psd.count <- non.psd.count + 1
    }
  }
  
  cat("non-PSD sample count:", non.psd.count, "\n")
  
  if (b < B){
    stop(paste0("Could not generate ", B, " PSD samples within max.iter = ", max.iter, "."))
  }
  
  if(output=="full"){
    return(Rhat)
  } else{
    return(t(rhat.keep))
  }
}

# Generate Correlation from bootstrapping data then compute correlation
boot_cor <- function(X, asy.n, B) {
  X <- as.matrix(X)
  p <- ncol(X)
  
  arr <- array(NA_real_, dim = c(p, p, B))
  for (b in seq_len(B)) {
    boot_id <- sample.int(asy.n, size = asy.n, replace = TRUE)
    arr[, , b] <- cor(X[boot_id, , drop = FALSE])
  }
  arr
}

# MCMC Cloud
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
    Psi_Rcurrent<- calc_psi(Rcurrent, n = asy.n)
    
    Rs<- rwaldcloud(Rcurrent, Psi= scale_factor*Psi_Rcurrent, B = 1) %>% drop()
    rs<- Rs[lower.tri(Rs)]
    
    Psi_Rs<- calc_psi(Rs, n= asy.n)
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
