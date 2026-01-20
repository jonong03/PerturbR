rm(list=ls());gc()
set.seed(123)
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm,rethinking, future.apply, parallel, parallelly)
wald.test <- function(Rpop, Rsample, alpha = 0.05, asy.n = 100000, fisherz = FALSE) {
  # Rpop     : Target (true) correlation matrix
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
    cov <- metaSEM::asyCov(R, n= asy.n)
    
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
    
    Psi <- metaSEM::asyCov(Rpop, n = asy.n)  # asymptotic covariance in r-space
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
  if(qr(Psi0)$rank != ps) {
    return(stop("Error: Covariance of R0 is not PSD"))
    
  }
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

d=10; n= 1e3
R1<- rlkjcorr(1,K=d,eta=2)
Rhat<- rlkjcorr(1,K=d,eta=2)
wald.test(R1, Rhat, fisherz= TRUE)

B=1000
Rb<- rwaldcloud(R1, n= 1e5, B=B, output="matrix")
out<- sapply(1:B, function(b) wald.test(R1, Rb[,,b], asy.n=1e5)) %>% t()
mean(out[,5])

# Plot density
p=d*(d-1)/2
x= out[,1]
hist(x, breaks = 30, probability = TRUE, main = "Histogram with Chi-square curve",xlab = "Wald-Statistic")
curve(dchisq(x, df= p),  add = TRUE, col = "red", lwd = 2)




Rb1<- lapply(1:B, function(i) cor(mvtnorm::rmvnorm(1e5, sigma = R1)))
out<- sapply(1:B, function(b) wald.test(R1, Rb1[[b]], asy.n=1e5)) %>% t()
mean(out[,5])
hist(out[,1])





R1<- rlkjcorr(1,K=d,eta=0.2)
Rb<- runifcloud(R1, n= n, B=B, output="matrix")
out<- sapply(1:dim(Rb)[3], function(b) wald.test(R1, Rb[,,b], asy.n=n)) %>% t()
mean(out[,5])


# plotly
R_to_coords <- function(R) {
  data.frame(
    r12 = R[1,2],
    r13 = R[1,3],
    r23 = R[2,3]
  )
}

d=3
R1<- diag(d)
R2<- rlkjcorr(1, K=d, eta= 8)
R3<- rlkjcorr(1, K=d, eta= 1)
R4<- rlkjcorr(1, K=d, eta= 0.2)

plot.corspace<- function(R0, B=1000, alpha=0.15, asy.n=1e5, main=""){
  r0_vech<- R_to_coords(R0)

  # Generate random samples
  rhat <- rwaldcloud(R0, B=B, n=asy.n, output= "matrix")
  rhat_vech <- t(sapply(1:B, function(b) R_to_coords(rhat[,,b])))
  rhat_vech <- data.frame(rhat_vech)
  colnames(rhat_vech)<- c("r12", "r13", "r23")
  
  # Perform Wald-test
  wald.out<- sapply(1:B, function(b)
    wald.test(Rpop= R0, Rsample = rhat[,,b], alpha = alpha, asy.n= asy.n)
    )
  outside.ellipsoid <- wald.out[5,]
  cov= 1-mean(outside.ellipsoid)

  inside = rhat_vech[which(outside.ellipsoid==0),]
  outside = rhat_vech[which(outside.ellipsoid==1),]
  
  N <- 3000
  samp1 <- rlkjcorr(N, K = 3, eta = 1) #eta parameter represents density of sampling
  df1 <- data.frame(r12=samp1[, 2, 1], r13=samp1[, 3, 1], r23=samp1[, 3, 2])
  
  title=main
  subtitle=paste0("Actual Coverage:", cov)
  fulltext= paste0(title,"<br><sup>", subtitle,"</sup>")
    
  pp<- plot_ly(
    data = inside,
    x = ~r12, y = ~r13, z = ~r23,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 4),
    name = "Inside samples"
  ) %>%
    add_trace(
      data = outside,
      x = ~r12, y = ~r13, z = ~r23,
      type = "scatter3d",
      mode = "markers",
      color = "black",
      marker = list(size = 4),
      name = "Outside samples"
    ) %>%
    add_trace(
      data = r0_vech, 
      x = ~r12, y = ~r13, z = ~r23,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 8, color = "red"),
      name = "R0"
    ) %>%
    add_trace(
      data = df1,
      x = ~r12, y = ~r13, z = ~r23,
      type = "scatter3d",
      mode = "markers",
      color = "lightgray",
      marker = list(size = 3),
      name = "Full space",
      visible = "legendonly"
    ) %>%
    layout(
      title = list(
        text = fulltext
        )
    ) 
    #layout(title = paste0("actual coverage:", cov))
  return(pp)
}
p1<- plot.corspace(R1, B=1000, alpha= 0.05, asy.n= 2e5, main="Identity Matrix")
p2<- plot.corspace(R2, B=1000, alpha= 0.05, asy.n= 2e5, main="eta= 8")
p3<- plot.corspace(R3, B=1000, alpha= 0.05, asy.n= 2e5, main="eta= 1")
p4<- plot.corspace(R4, B=1000, alpha= 0.05, asy.n= 2e5, main="eta= 0.2")

pacman::p_load(htmlwidgets)
saveWidget(p1, "docs/plots/p1.html", selfcontained = FALSE)
saveWidget(p2, "docs/plots/p2.html", selfcontained = FALSE)
saveWidget(p3, "docs/plots/p3.html", selfcontained = FALSE)
saveWidget(p4, "docs/plots/p4.html", selfcontained = FALSE)


plot.corspace.unif<- function(R0, B=1000, alpha=0.15, asy.n=1e5, main=""){
  r0_vech<- R_to_coords(R0)
  
  # Generate random samples
  rhat <- runifcloud(R0, B=B, n=asy.n, alpha=alpha, output= "matrix")
  rhat_vech <- t(sapply(1:B, function(b) R_to_coords(rhat[,,b])))
  rhat_vech <- data.frame(rhat_vech)
  colnames(rhat_vech)<- c("r12", "r13", "r23")
  
  # Perform Wald-test
  wald.out<- sapply(1:B, function(b)
    wald.test(Rpop= R0, Rsample = rhat[,,b], alpha = alpha, asy.n= asy.n)
  )
  outside.ellipsoid <- wald.out[5,]
  cov= 1-mean(outside.ellipsoid)
  
  inside = rhat_vech[which(outside.ellipsoid==0),]
  outside = rhat_vech[which(outside.ellipsoid==1),]
  
  N <- 3000
  samp1 <- rlkjcorr(N, K = 3, eta = 1) #eta parameter represents density of sampling
  df1 <- data.frame(r12=samp1[, 2, 1], r13=samp1[, 3, 1], r23=samp1[, 3, 2])
  
  title=main
  subtitle=paste0("Actual Coverage:", cov)
  fulltext= paste0(title,"<br><sup>", subtitle,"</sup>")
  
  pp<- plot_ly(
    data = inside,
    x = ~r12, y = ~r13, z = ~r23,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 4),
    name = "Inside samples"
  ) %>%
    add_trace(
      data = outside,
      x = ~r12, y = ~r13, z = ~r23,
      type = "scatter3d",
      mode = "markers",
      color = "black",
      marker = list(size = 4),
      name = "Outside samples"
    ) %>%
    add_trace(
      data = r0_vech, 
      x = ~r12, y = ~r13, z = ~r23,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 8, color = "red"),
      name = "R0"
    ) %>%
    add_trace(
      data = df1,
      x = ~r12, y = ~r13, z = ~r23,
      type = "scatter3d",
      mode = "markers",
      color = "lightgray",
      marker = list(size = 3),
      name = "Full space",
      visible = "legendonly"
    ) %>%
    layout(
      title = list(
        text = fulltext
      )
    ) 
  #layout(title = paste0("actual coverage:", cov))
  return(pp)
}
p1unif<- plot.corspace.unif(R1, B=1000, alpha= 0.05, asy.n= 2e5, main="Identity Matrix")
p2unif<- plot.corspace.unif(R2, B=1000, alpha= 0.05, asy.n= 2e5, main="eta= 8")
p3unif<- plot.corspace.unif(R3, B=1000, alpha= 0.05, asy.n= 2e5, main="eta= 1")
p4unif<- plot.corspace.unif(R4, B=1000, alpha= 0.05, asy.n= 2e5, main="eta= 0.2")

saveWidget(p1unif, "docs/plots/p1unif.html", selfcontained = FALSE)
saveWidget(p2unif, "docs/plots/p2unif.html", selfcontained = FALSE)
saveWidget(p3unif, "docs/plots/p3unif.html", selfcontained = FALSE)
saveWidget(p4unif, "docs/plots/p4unif.html", selfcontained = FALSE)

# Simulation
# ---- build parameter grid ----

####
eta_vals   <- c(0.1, 1, 8)
alpha_vals <- c(0.05, 0.10)
asy_n_vals <- c(1e6)
iter_vals  <- 300
nG         <- 50             # number of graphs / replicates of graphs
fisherz    <- c(FALSE, TRUE)
dimension  <- c(3, 5, 10)

grid <- expand.grid(
  eta   = eta_vals,
  alpha = alpha_vals,
  asy.n = asy_n_vals,
  iter  = iter_vals,
  nG     = seq_len(nG),
  fisherz= fisherz, 
  dimension = dimension,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$nParam = grid$dimension * (grid$dimension -1)/2
#grid$alpha = grid$alpha/grid$nParam
base_seed <- 123
grid$seed <- base_seed + seq_len(nrow(grid))
grid$Rtarget <- lapply(seq_len(nrow(grid)), function(i) {
  set.seed(grid$seed[i])
  rlkjcorr(1, grid$dimension[i], eta = grid$eta[i])
})
simulate_coverage_joint <- function(R0, alpha = 0.05, asy.n = 10000, iter = 100, fisherz = FALSE) {
  
  Rb <- tryCatch(
    rwaldcloud(R0, n = asy.n, B = iter, output = "matrix"),
    error = function(e) NA
  )
  
  # If rwaldcloud failed, return NA
  if (all(is.na(Rb))) return(NA_real_)
  
  validB <- dim(Rb)[3]
  out <- sapply(1:validB, function(b) { wald.test(R0, Rb[, , b], asy.n = asy.n, alpha= alpha)}) %>% t()
  
  cov <- 1 - mean(out[, 5], na.rm = TRUE)
  return(cov)
}



{
  plan(multisession, workers = max(1, parallelly::availableCores()- 1))
  res_joint <- future_sapply(seq_len(nrow(grid)), function(i) {
    res <- simulate_coverage_joint(
      R0 = grid$Rtarget[[i]],
      alpha   = grid$alpha[i],
      asy.n   = grid$asy.n[i],
      iter    = grid$iter[i],
      fisherz = as.logical(grid$fisherz[i])
    )
    #res.out= rowMeans(res, na.rm=TRUE)
    return(res) 
  }, future.seed = TRUE)
  grid.joint<- as.data.table(grid)
  grid.joint<- cbind(grid.joint, coverage= res_joint)
  grid.joint[,Rtarget:= NULL]
  
  plan(sequential)
}
grid.joint[eta==0.1 & alpha==0.05,]

out<- grid.joint[,{
  qs <- quantile(coverage, c(0, .25, .5, .75, 1), na.rm = TRUE)
  .(
    Mean = sprintf("%.3f", mean(coverage, na.rm = TRUE)),
    #Min  = sprintf("%.3f", qs[1]),
    Q1   = sprintf("%.3f", qs[2]),
    #Q2   = sprintf("%.3f", qs[3]),
    Q3   = sprintf("%.3f", qs[4]),
    countNA = sum(is.na(coverage))
    #Max  = sprintf("%.3f", qs[5])
  )
}, by=.(eta, dimension, nParam, fisherz, asy.n, alpha)] 

saveRDS(out, "docs/objects/sim1.rds")



## CDA Sim
getwd()
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R")
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/EST/CDAfunctions.R")

alpha1=0.1
n0<- 1e6
d<- 10
iter=500

R0 <- rlkjcorr(1, K=d, eta=3)
fit_obj<- pcalg::pc(suffStat= list(C = R0, n = n0), indepTest = gaussCItest, p= d, alpha = as.numeric(0.05) )
G0<- true_adj<- t(as(fit_obj, "matrix"))
sum(G0)

Rb<- runifcloud(R0, n= n0, B=iter, alpha= alpha)
perf_Gb<-sapply(1:iter, function(b) est_pc(R=Rb[,,b], n=n0, alpha=alpha)) %>% t()
apply(perf_Gb, 2, summary) %>% round(.,1)

Gb<- lapply(1:iter, function(b){
  fit_obj<- pcalg::pc(suffStat= list(C = Rb[,,b], n = n0), indepTest = gaussCItest, p= d, alpha = as.numeric(0.05) )
  return( t(as(fit_obj, "matrix")) )
})


simulate_pc_performance <- function(R0, alpha, n0, iter, fisherz = FALSE, alpha2= 0.05) {
  
  d = nrow(R0)
  # Baseline PC graph from R0 (treat R0 as population correlation)
  fit0<- pcalg::pc(suffStat= list(C = R0, n = n0), indepTest = gaussCItest, p= d, alpha = alpha2 )
  G0<- true_adj<- t(as(fit0, "matrix"))
  edges_G0 <- sum(G0)
  
  # Generate cloud of correlations around R0
  Rb <- tryCatch(
    runifcloud(R0, n = n0, B = iter, alpha = alpha),
    error = function(e) NA
  )
  
  validB <- dim(Rb)[3]
  
  # Evaluate PC per cloud draw
  perf_list <- sapply(seq_len(validB), function(b) {
    out <- tryCatch(
      est_pc(R = Rb[, , b], n = n0, alpha = alpha2),
      error = function(e) NA
    )
    out
  })
  
  return(t(perf_list))
}

R0 <- rlkjcorr(1, K=12, eta=0.6)
out1<- simulate_pc_performance(R0, alpha= 0.05, n0= n0, iter=5)
apply(out1,2,summary)

out2<- simulate_pc_performance(R0, alpha= 0.01, n0= n0, iter=5)
apply(out2,2,summary)

saveRDS(list(out1 = out1, out2 = out2), "docs/objects/pc_perf.rds")


## Scale Up Sim

####
eta_vals   <- c(0.1, 1, 8)
alpha_vals <- c(0.05, 0.10)
asy_n_vals <- c(1e6)
iter_vals  <- 300
nG         <- 50             # number of graphs / replicates of graphs
fisherz    <- c(FALSE, TRUE)
dimension  <- c(5, 10)

grid <- expand.grid(
  eta   = eta_vals,
  alpha = alpha_vals,
  asy.n = asy_n_vals,
  iter  = iter_vals,
  nG     = seq_len(nG),
  fisherz= fisherz, 
  dimension = dimension,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$nParam = grid$dimension * (grid$dimension -1)/2
#grid$alpha = grid$alpha/grid$nParam
base_seed <- 123
grid$seed <- base_seed + seq_len(nrow(grid))
grid$Rtarget <- lapply(seq_len(nrow(grid)), function(i) {
  set.seed(grid$seed[i])
  rlkjcorr(1, grid$dimension[i], eta = grid$eta[i])
})

