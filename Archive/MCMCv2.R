# New: MCMC-Style Sampling to Obtain Confidence Set of R

# Step 0: Initialize chain R with a starting value: S (observed R)
# Step 1: Draw a proposal R* from wald cloud of R^(l-1) : H ~ W(R^(l-1))
# Step 2: Decide to accept or reject R* based on a state-dependent MH ratio 
# Step 3: If accept, R^l <- R* 
#         If reject, R^l <- R^(l-1) 
# Repeat Step 1,2,3

# Check: acceptance rate, calibration (p-value)

rm(list=ls()); gc()
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr, plotly)

use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
py_config() 
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)

# Acceptance Rate: ~55% (p=3), ~ 22% (p=5), ~2% (p=8). Sample Size = 500
# ~0.01% (p=10), sample size = 50K


# Prep: Generate trueR
p=3L; ad=2L;  NCHAIN= 10000L; N= 50000L
df<- m<- p*(p-1)/2  # degree of freedom
Target <- er_dag_py(p = p, ad = ad, n = N, K = 1L, seed = 144L)
G0 <- Target$G
R0 <- Target$R ; kappa(R0)

asy.n= 5000L
sampleid<- sample(N, asy.n)
X <- Target$X[1,sampleid,]
S <- cor(X); kappa(S)


# S = rlkjcorr(1, K=3, eta = 5)
# Step 0: Initiation of Chain
RCHAIN<- array(NA, dim= c(p,p,NCHAIN))
RCHAIN[,,1] <- S 

type = 1
scale_factor = 0.4
accept <- 0

for (i in 2: NCHAIN){
  Rcurrent<- RCHAIN[,,i-1]
  rcurrent<- Rcurrent[lower.tri(Rcurrent)]
  
  # Step 1: Draw R*
  
  Rs<- rwaldcloud(Rcurrent, n = asy.n, B = 1, scale= scale_factor) %>% drop()
  rs<- Rs[lower.tri(Rs)]
  
  # Step 2: Acceptance Ratio
  Psi_Rcurrent<- metaSEM::asyCov(Rcurrent, n= asy.n)
  Psi_Rs<- metaSEM::asyCov(Rs, n= asy.n)
  
  partA<- dmvnorm(x= rcurrent, mean = rs, sigma= Psi_Rs, log= T)
  partB<- dmvnorm(x= rs, mean= rcurrent, sigma= Psi_Rcurrent, log=T)
  # Chi-square/ F-distribution
  partC<- wald.test(Rpop=Rs, Rsample= S, alpha= 0.05, asy.n= asy.n)  # Rstar is at the center
  #cat(partA,"\n")
  IR   <- partC$pval >0.05
  # partC$pval < 0.05  # if p-val< 0.05, reject, outside of 95% CR
  logMH = partA - partB 
  
  ###NEW
  if(type==2){
    D1 <- t(rcurrent - rs) %*% solve(Psi_Rs) %*% (rcurrent-rs)
    D2 <- t(rcurrent - rs) %*% solve(Psi_Rcurrent) %*% (rcurrent-rs)
    detD1 <- determinant(Psi_Rs, logarithm = T)$modulus
    detD2 <- determinant(Psi_Rcurrent, logarithm = T)$modulus
    logMH = (-0.5* (D1-D2) - 0.5*(detD1 - detD2))[1]
  }
  
  # Step 3: accept or reject
  if(IR==TRUE && logMH > log(runif(1))){ # accept Rs
    RCHAIN[,,i]<- Rs
    accept <- accept + 1
  } else {
    RCHAIN[,,i]<- Rcurrent
  }
  cat("Accept:", accept/i," ")
}

acceptrate= accept/NCHAIN
acceptrate


burnin<- 1:(NCHAIN*0.20)

dim(RCHAIN)
Qi <- sapply((NCHAIN*0.20+1):NCHAIN, function(i){
  rchain<- RCHAIN[,,i]
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
summary(Qi)
hist(Qi)

qchisq(p=0.05, d=3, lower.tail = F)


####
eigenvec1<- eigen(Psi_Rs)$vectors
eigenval1<- diag(ncol(Psi_Rs))
diag(eigenval1) <- eigen(Psi_Rs)$values
eigenvec1 %*% eigenval1 %*% t(eigenvec1)
solve(Psi_Rs)

L1<- eigenvec1 %*% sqrt(eigenval1)
# L1 %*% t(L1)

y1<- solve(L1) %*% (rcurrent- rs)
qchisq(p = 0.05, df=3)
####


# Check Reverse Containment: Each R should be compatible with S
out<- sapply(c(1:NCHAIN)[-burnin], function(i) wald.test(Rsample=S, Rpop= RCHAIN[,,i], alpha= 0.05, asy.n= asy.n)) %>% t()
pval<- sapply(1:nrow(out), function(i) out[,4][[i]])
hist(pval)
rcr<- mean(pval > 0.05)
rcr
summary(pval)

# Generate wald cloud, inverse test (reverse containment)
waldCR <- rwaldcloud(R0=S, n = asy.n, B=NCHAIN)
pval_wald<- sapply(1:NCHAIN, function(i) {
  out<- wald.test(Rsample=S, Rpop= waldCR[,,i], alpha= 0.05, asy.n= asy.n)
  out$pval
  })
hist(pval_wald)
waldCR_in<- waldCR[,,which(pval_wald > 0.05)]
pval_wald2<- sapply(1:sum(pval_wald> 0.05), function(i) {
  out<- wald.test(Rsample=S, Rpop= waldCR_in[,,i], alpha= 0.05, asy.n= asy.n)
  out$pval
})
mean(pval_wald > 0.05)
mean(pval_wald2 > 0.05)
rcr_wald<- mean(pval_wald2 > 0.05)
rcr_wald

fullspace<- rlkjcorr(n=NCHAIN*100, K=3)
dim(fullspace)
pval_fullspace<- sapply(1:dim(fullspace)[1], function(i){
  R<- fullspace[i,,]
  out<- wald.test(Rsample= S, Rpop= R, alpha = 0.05, asy.n= asy.n)
  out[4]
})
pval_fullspace<- sapply(pval_fullspace, function(i) i)
fullspace_in <- fullspace[which(pval_fullspace> 0.05),,]
dim(fullspace_in)

hist(pval_fullspace)
sum(pval_fullspace> 0.05, na.rm=T)
8/40000



# Plot --------------------------------------------------------------------
{
  sampleNCHAIN= sample(NCHAIN*100, NCHAIN)
  fullspace_df<- data.frame(
    r12 = fullspace[sampleNCHAIN,1,2],
    r13 = fullspace[sampleNCHAIN,1,3],
    r23 = fullspace[sampleNCHAIN,2,3]
  )
  fullspace_in_df<- data.frame(
    r12 = fullspace_in[,1,2],
    r13 = fullspace_in[,1,3],
    r23 = fullspace_in[,2,3]
    
  )
  chain_df <- data.frame(
    r12 = RCHAIN[1,2,-burnin],
    r13 = RCHAIN[1,3,-burnin],
    r23 = RCHAIN[2,3,-burnin]
  )
  S_point <- data.frame(
    r12 = S[1,2],
    r13 = S[1,3],
    r23 = S[2,3]
  )
  R_point <- data.frame(
    r12 = R0[1,2],
    r13 = R0[1,3],
    r23 = R0[2,3]
  )
  waldCR_df<- data.frame(
    r12= waldCR_in[1,2,],
    r13= waldCR_in[1,3,],
    r23= waldCR_in[2,3,]
  )
  
  
resultplot<-plot_ly() %>%
    add_markers(
      data = fullspace_df,
      x = ~r12,
      y = ~r13,
      z = ~r23,
      marker = list(size = 2, color = "gray"),
      opacity = 0.2,
      name = "Full space",
      visible = "legendonly"
      ) %>%
  add_markers(
    x = c(0,0),
    y = c(0,0),
    z = c(0,0),
    marker = list(size = 8, color ="black"),
    name = "Identity Matrix"
  ) %>%
  add_markers(
    data = R_point,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    marker = list(size = 8, color = "orange"),
    name = "R (unknown)"
  ) %>%
  add_markers(
    data = S_point,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    marker = list(size = 8, color = "red"),
    name = "S (observed)"
  ) %>%
  add_markers(
    data = waldCR_df,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    marker = list(size = 2, color = "darkgreen"),
    opacity = 0.4,
    name = "Wald space of S - Rev.Containment"
  )%>%
    add_markers(
      data = chain_df,
      x = ~r12,
      y = ~r13,
      z = ~r23,
      marker = list(size = 2, color = "blue"),
      opacity = 0.4,
      name = "MCMC samples"
    ) %>%
  add_markers(
    data = fullspace_in_df,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    marker = list(size = 3, color = "gold"),
    name = "Full space of 3x3 - Rev.Containment"
  )%>%
    layout(
      scene = list(
        xaxis = list(title = "r12"),
        yaxis = list(title = "r13"),
        zaxis = list(title = "r23")
      ),
      legend = list(
        itemsizing = "constant"
      ),
      title = list(
        text = paste(
          "MCMC Cloud vs Wald Cloud (n=",asy.n,")",
          paste0("<br><sup> Reverse Containment Rate- MCMC:", round(rcr*100,1),"%; Wald Cloud: ", round(rcr_wald*100,1), "%<br> MCMC Acceptance Rate:", round(acceptrate,2), "</sup>")
        ),
        x = 0.5, y= 0.95   # center the title
      )
    )
}

resultplot


library(htmlwidgets)

savedir<- "~/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/docs/plots/ReverseContainment_MCMC/"
saveWidget(resultplot, paste0(savedir,"mcmc_cloud_n",asy.n,".html"), selfcontained = TRUE)

asy.n

# Diagnostic --------------------------------------------------------------
png(paste0(savedir,"trace_plots_n",asy.n,".png"), width = 800, height = 800)
par(mfrow=c(3,1))
plot(RCHAIN[1,2,-burnin], type = "l")
plot(RCHAIN[1,3,-burnin], type = "l")
plot(RCHAIN[2,3,-burnin], type = "l")
dev.off()

# ACF
png(paste0(savedir,"acf_plots_n",asy.n,".png"), width = 950, height = 400)
par(mfrow=c(1,3))
idx <- seq(1,NCHAIN, by=50)
acf(RCHAIN[1,2,-burnin], lag.max = 50)
acf(RCHAIN[1,3,-burnin], lag.max = 50)
acf(RCHAIN[2,3,-burnin], lag.max = 50)
dev.off()

acf(RCHAIN[1,2,idx], lag.max = 50)
acf(RCHAIN[1,3,idx], lag.max = 50)
acf(RCHAIN[2,3,idx], lag.max = 50)


# marginal distribution

par(mfrow=c(2,3))
hist(RCHAIN[1,2,-burnin])
hist(RCHAIN[1,3,-burnin])
hist(RCHAIN[2,3,-burnin])

hist(waldCR_df[,1])
hist(waldCR_df[,2])
hist(waldCR_df[,3])

