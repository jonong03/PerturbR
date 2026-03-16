# MCMC-Style Sampling to Obtain Confidence Set of R

# Step 0: Initialize chain with a starting value: S (observed R)
# Step 1: Draw a proposal H from wald cloud of S : H ~ W(S)
# Step 2: Test H0: S= H (centered at H), reject when p-val < 0.05
# Step 3: If accept, add H to the chain
#         If reject, no updates to the chain
# Step 4: Draw another proposal H* ~ W(S*), where S* is the latest value inthe chain, and repeat the remaining step 2 and 3.


rm(list=ls()); gc()
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr)

use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
py_config() 
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)

# 1: Generate trueR
p=3L; ad=2L; asy.n= 500L; B= 5000L
df<- m<- p*(p-1)/2  # degree of freedom
Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 1L)
G0 <- Target$G
R0 <- Target$R
X <- Target$X[1,,]

S <- cor(X)
S; R0

NCHAIN = B
RHbag<- array(NA, dim= c(p,p,NCHAIN))
accept <- 0

# Initialize Chain
RHbag[,,1]<- S

for(i in 2:NCHAIN){
  
  H <- rwaldcloud(R0=RHbag[,,i-1], asy.n)[,,1]  # propose candidate H based on wald cloud of previous RH
  #covH check
  pval <- wald.test(Rpop= H, Rsample= S, asy.n=asy.n, alpha = 0.05)[4]
  
  if(pval > 0.05){
    RHbag[,,i]<- H 
    accept <- accept + 1
  } else {
    RHbag[,,i]<- RHbag[,,i-1]
  }
} ++

dim(RHbag)
accept/NCHAIN

out1<- sapply(1:NCHAIN, function(i) wald.test(Rpop= RHbag[,,i], Rsample= S, asy.n = asy.n)) %>% t()
head(out1)
summary(out1[,4])

# Take unique set?


out2<- sapply(1:NCHAIN, function(i) wald.test(Rpop= S, Rsample= RHbag[,,i], asy.n = asy.n)) %>% t()
head(out2)
summary(out2[,4])
mean(out2[,4]<0.05)

kappa(S)

# How does this compared with forward 95% CI approach, 
S2 <- rwaldcloud(R0=S, n= asy.n, B=B)
validCI<- sapply(1:B, function(i) wald.test(Rpop = S, Rsample = S2[,,i], alpha = 0.05, asy.n = asy.n))
mean(validCI[4,]<0.05)
inside<- which(validCI[4,]>0.05)
S2inside<- S2[,,inside]
dim(S2inside)

# Put into a data frame
df2 <- data.frame(
  r12 = S2inside[1,2,],
  r13 = S2inside[1,3,],
  r23 = S2inside[2,3,]
)





library(plotly)

# Suppose your array is called R_array
# dim(R_array) should be c(3, 3, 1000)

# Extract the off-diagonal correlations
r12 <- RHbag[1, 2, ]
r13 <- RHbag[1, 3, ]
r23 <- RHbag[2, 3, ]

# Put into a data frame
df <- data.frame(
  r12 = r12,
  r13 = r13,
  r23 = r23
)
S_point <- data.frame(
  r12 = S[1,2],
  r13 = S[1,3],
  r23 = S[2,3]
)
# Basic 3D scatter plot
plot_ly() %>%
  # Cloud of sampled correlation matrices
  add_trace(
    data = df,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 3,
      color = "grey",
      opacity = 0.6
    ),
    name = "R_H"
  ) %>%
  
  # Add S as red dot
  add_trace(
    data = S_point,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 8,
      color = "red"
    ),
    name = "S"
  ) %>%
  # Add S as red dot
  add_trace(
    x = R0[1,2],
    y = R0[1,3],
    z = R0[2,3],
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 8,
      color = "green"
    ),
    name = "R"
  ) %>%
  layout(
    title = "3D Correlation Space",
    scene = list(
      xaxis = list(title = "r12", range = c(-1,1)),
      yaxis = list(title = "r13", range = c(-1,1)),
      zaxis = list(title = "r23", range = c(-1,1))
    )
  ) %>%
  add_trace(
    data = df2,
    x = ~r12,
    y = ~r13,
    z = ~r23,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 3,
      color = "lightblue",
      opacity = 0.6
    ),
    name = "95%CI(S)"
  )


wald.with.R_H<- sapply(1:dim(RHbag)[3], function(i) wald.test(Rpop= S, Rsample = RHbag[,,i], alpha= 0.05, asy.n= asy.n)) %>% t()
wald.with.CI_S<- sapply(1:dim(S2inside)[3], function(i) wald.test(Rpop= R0, Rsample = S2inside[,,i], alpha= 0.05, asy.n= asy.n)) %>% t()
i=1
wald.with.R_H<- sapply(1:dim(RHbag)[3], function(i) wald.test(Rsample= S, Rpop = RHbag[,,i], alpha= 0.05, asy.n= asy.n)) %>% t()
wald.with.CI_S<- sapply(1:dim(S2inside)[3], function(i) wald.test(Rsample= S, Rpop = S2inside[,,i], alpha= 0.05, asy.n= asy.n)) %>% t()

mean(wald.with.R_H[,4]< 0.05)
mean(wald.with.CI_S[,4]< 0.05)
hist(wald.with.R_H[,4])



pacman::p_load(mvtnorm)
s<- S[upper.tri(S)]
h<- H[upper.tri(H)]
sigH<- metaSEM::asyCov(H, n= asy.n)
sigS<- metaSEM::asyCov(S, n= asy.n)
dmvnorm(h, s, sigS)
dmvnorm(s, h, sigH)
