rm(list=ls());gc()
source("docs/WorkingCode/Functions.R")

pacman::p_load(metaSEM, dplyr, data.table, mvtnorm,rethinking, future.apply, parallel, parallelly)

set.seed(123)
d=10; n= 1e3
R1<- rlkjcorr(1,K=d,eta=2)
Rhat<- rlkjcorr(1,K=d,eta=2)
wald.test(R1, Rhat, fisherz= TRUE)

B=1000
Rb<- rwaldcloud(R1, n= 1e5, B=B, output="matrix")
out<- sapply(1:B, function(b) wald.test(R1, Rb[,,b], asy.n=1e5)) %>% t()
mean(out[,5])

# Check theoretical and actual density

R_d3<- rlkjcorr(1,K=3, eta=1)
R_d5<- rlkjcorr(1, K=5, eta=1)
R_d10<- rlkjcorr(1,K=10, eta=1)
R_d20<- rlkjcorr(1, K=20, eta=1)

## Simulate actual and theoretical density of test statistic
density.simulate<- function(R0, B=1000, asy.n=100000){
  d= ncol(R0)
  p=d*(d-1)/2
  Rb<- rwaldcloud(R0, n= asy.n, B=B, output="matrix")
  out<- sapply(1:B, function(b) wald.test(Rpop= R0, Rsample= Rb[,,b], asy.n= asy.n))
  x<- out[1,] # test statistic
  
  hist(x, breaks = 30, probability = TRUE, main = paste0("#Nodes (d):",d, "; degree-of-freedom (p)=",p) ,xlab = "Wald-Statistic")
  curve(dchisq(x, df= p),  add = TRUE, col = "red", lwd = 2)
}

par(mfrow=c(2,2), oma = c(0, 0, 3, 0))

density.simulate(R_d3)
density.simulate(R_d5)
density.simulate(R_d10)
density.simulate(R_d20)

mtext("Empirical histogram with Chi-square curve (red)", outer = TRUE, cex = 1.4, line = 1)







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

