# Sim5_ Reverse Containment

# 1: Generate trueR
# 2: Make perturbations around true R: Many Rhats
# 3: Check if Rhats contain the trueR


# Sim1:
# 1a: When variance of rho0 is known, use Psi(rho0)
# 1b: When variance of rho0 is unknown, use sample estimate of Psi(rho0) (Sample covariance matrix).
# 1c: When both mean and variance of rho0 is unknown, 
#     use sample estimate of both rho0: rho_hat, and sampel estimate of Psi(rho0) (Sample covariance matrix). 
#     Do we call this the "empirical distribution" of mahalanobis distance.
# Distribution of resuling mahalanobis distance in each cases above:
# 1a: Original definition of mahalanobis distance: Chi-square with df m (sum of squares of m independent normally distributed variables: rho).  
# 1b: 


rm(list=ls()); gc()
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr)


use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
py_config() 
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)


######## Test #######
# 1: Generate trueR
p=30L; ad=2L; asy.n= 50000L; B= 1000
df<- m<- p*(p-1)/2  # degree of freedom
Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 1L)
G0 <- Target$G
R0 <- Target$R

p=30L; ad=2L; asy.n= 50000L; B= 1000
df<- m<- p*(p-1)/2  # degree of freedom

R0 <- rlkjcorr(n=1, K= p, eta = 1)  # Prefer this, and test over different eta values, make sure R0 is well-defined
Psi_R0<- metaSEM::asyCov(R0, n= asy.n)
qr(Psi_R0)$rank==ncol(Psi_R0)  # TRUE = well-defined


# 2: Generate Perturbations around R0
Rhat <- rwaldcloud(R0, n= asy.n, B=B)   # returns array
Psi_Rhat <- lapply(1:B, function(i) metaSEM::asyCov(x= Rhat[,,i], n= asy.n))

## Check which Psi_Rhat is ill-defined
Psi_Rhat_ill <- sapply(Psi_Rhat, function(Q) qr(Q*asy.n)$rank != ncol(Q))
Psi_Rhat_ill %>% sum


# 3: Check if Rhats contain the trueR using mahalanobis distance
# 3-1: Compute theoretical variance for each trueR

waldstat<- sapply(1:B, function(i){
  ps<- ncol(Rpop)
  Rsample <- Rhat[,,i]
  
  Rpop_vec<- R0[lower.tri(R0)]
  Rsample_vec <- Rsample[lower.tri(Rsample)]
  Psi_Rhat.i<- metaSEM::asyCov(x= Rsample, n= asy.n, cor.analysis= T)
  
  ## --- Mahalanobis distance ---
  distance = mahalanobis(x= Rsample_vec,center = Rpop_vec,cov = Psi_Rhat[[i]])
  return(distance)
})

hist(waldstat, probability = TRUE, col = "lightgray", border = "white", main = "Histogram with Theoretical Density", breaks= 20)


alpha= 0.05
cricval<- qchisq(1 - alpha, df = df)
pval<- as.numeric(pchisq(waldstat, df= df, lower.tail = F))
reject= as.logical(waldstat >= cricval)
out = cbind(T=waldstat, df= df, cricval= as.numeric(cricval), pval= pval, reject= reject) 
mean(reject)


par(mfrow=c(2,2))
hist(waldstat, main= "Distribution of wald distances")
qqplot(
  qchisq(ppoints(length(waldstat)), df = m),
  sort(waldstat),
  xlab = "Theoretical Chi-square Quantiles",
  ylab = "Observed Quantiles",
  main = "QQ plot of Mahalanobis distances"
)
abline(0,1,col="red",lwd=2)

hist(pval, breaks = 20, main= "Distribution of p-value")

qqplot(
  qunif(ppoints(length(pval))),  # theoretical uniform quantiles
  sort(pval),                    # observed quantiles
  xlab = "Theoretical Quantiles (Uniform)",
  ylab = "Observed Quantiles",
  main = "QQ plot of p-values"
)
abline(0, 1, col = "red", lwd = 2)





###############################

## settings
p <- 40L
B <- 1000L
eta <- 5  # 0.5, 5
n_vec <- c(500L, 5000L, 50000L)
df <- p * (p - 1) / 2
alpha <- 0.05

## same system: generate R0 only once
R0 <- rlkjcorr(n = 1, K = p, eta = eta)
R0_vec <- R0[lower.tri(R0)]
Psi_R0<- metaSEM::asyCov(R0, n= min(n_vec))
qr(Psi_R0)$rank==ncol(Psi_R0)  # TRUE = well-defined

## store results
res <- vector("list", length(n_vec))
names(res) <- n_vec

for (j in seq_along(n_vec)) {
  asy.n <- n_vec[j]
  
  ## generate perturbations around same R0
  Rhat <- rwaldcloud(R0, n = asy.n, B = B)
  
  Psi_Rhat <- lapply(1:B, function(i) metaSEM::asyCov(x= Rhat[,,i], n= asy.n))
  Psi_Rhat_ill <- sapply(Psi_Rhat, function(Q) qr(Q*asy.n)$rank != ncol(Q))

  ## Wald distances
  
  waldstat <- sapply(1:B, function(i) {
    Rsample <- Rhat[,,i]
    Rsample_vec <- Rsample[lower.tri(Rsample)]
    Psi_i <- Psi_Rhat[[i]]

    mahalanobis(
      x = Rsample_vec,
      center = R0_vec,
      cov = Psi_i
    )
  })
  
  ## p-values
  pval <- pchisq(waldstat, df = df, lower.tail = FALSE)
  
  res[[j]] <- list(
    n = asy.n,
    waldstat = waldstat,
    pval = pval,
    ill = Psi_Rhat_ill
  )
}

# Reverse Containment Rate for each setting
rate= sapply(1:3, function(i) mean(res[[i]]$pval<= alpha))
ill.defined.sample.rate = sapply(1:3, function(i) res[[i]]$ill %>% mean) # how many pertubed samples has non-full rank covariance?

## plotting
{
  layout(
    matrix(1:12, nrow=3, byrow=TRUE),
    heights=c(1.5,3,3)
  )
  par(mar = c(4, 4, 3, 1), oma=c(0,0,3,0))
  
  #par(mfrow = c(4, 4), mar = c(4, 4, 3, 1), oma=c(0,0,3,0))
  cols <- c("red", "blue", "darkgreen")

  ## row 1: Reverse containment rate
  for(i in 1:3){
    plot.new()
    text(
      0.5, 0.5,
      paste0("n = ", n_vec[i], "\nReverse containment%: ", round(1-rate[i],4),"\n non-PSD cov%: ", ill.defined.sample.rate[i]),
      cex = 1.6,
      font = 1
    )
  }
  plot.new()   # empty last panel
  
  
  ## -------------------
  ## row 3: Wald distance
  ## -------------------
  for (j in 1:3) {
    hist(
      res[[j]]$waldstat,
      probability = TRUE,
      breaks = 20,
      col = "lightgray",
      border = "white",
      main = paste("Wald distance\nn =", res[[j]]$n),
      xlab = "Wald distance"
    )
    curve(dchisq(x, df = df), add = TRUE, col = "red", lwd = 2)
  }
  
  ## overlapped QQ plot for Wald distance
  qq_x <- qchisq(ppoints(B), df = df)
  plot(
    qq_x, sort(res[[1]]$waldstat),
    type = "p",
    col = cols[1],
    pch = 16,
    cex = 0.6,
    xlab = "Theoretical chi-square quantiles",
    ylab = "Observed quantiles",
    main = "QQ plot of Wald distance"
  )
  points(qq_x, sort(res[[2]]$waldstat), col = cols[2], pch = 16, cex = 0.6)
  points(qq_x, sort(res[[3]]$waldstat), col = cols[3], pch = 16, cex = 0.6)
  legend("topleft", legend = paste("n =", n_vec), col = cols, pch = 16, bty = "n")
  
  ## -------------------
  ## row 2: p-values
  ## -------------------
  for (j in 1:3) {
    hist(
      res[[j]]$pval,
      probability = TRUE,
      breaks = 20,
      col = "lightgray",
      border = "white",
      main = paste("p-values\nn =", res[[j]]$n),
      xlab = "p-value"
    )
    abline(h = 1, col = "red", lwd = 2)
  }
  
  ## overlapped QQ plot for p-values
  qq_u <- ppoints(B)
  plot(
    qq_u, sort(res[[1]]$pval),
    type = "p",
    col = cols[1],
    pch = 16,
    cex = 0.6,
    xlab = "Theoretical uniform quantiles",
    ylab = "Observed quantiles",
    main = "QQ plot of p-values"
  )
  points(qq_u, sort(res[[2]]$pval), col = cols[2], pch = 16, cex = 0.6)
  points(qq_u, sort(res[[3]]$pval), col = cols[3], pch = 16, cex = 0.6)
  abline(0, 1, col = "black", lwd = 2, lty = 2)
  legend("topleft", legend = paste("n =", n_vec), col = cols, pch = 16, bty = "n")
  mtext(paste0(c("nodes:", "eta:","alpha:","Iter:"), c(p, eta, alpha,B), collapse=", "), outer=TRUE, cex=1.4)
  
}


###############################
# Repeat for the whole space

# 1- use rlkj to construct the whole space using eta =1 

# 2- Treat each point as truth (R0_i), construct perturbed samples for each R0_i: Rhat_ij

# 3- Check reverse containment rate: plot it as color. 
#    NA / gray, if a R0_i has non-full rank covariance (therefore unable to generate samples) 

set.seed(123)

## settings
p <- 40L
B <- 1000L
S <- 100L                  # number of different R0's
eta <- 1
n_vec <- c(500L, 5000L, 50000L)
df <- p * (p - 1) / 2
alpha <- 0.05
R0_all <- rlkjcorr(n = S, K = p, eta = eta)


## store results for all R0 systems
res_all <- vector("list", S)

#### Sequential Run
for (s in 1:S) {
  cat("Running R0 system", s, "of", S, "\n")
  
  ## generate one new R0
  R0 <- R0_all[s,,]
  R0_vec <- R0[lower.tri(R0)]
  
  ## store results across n for this R0
  res_n <- vector("list", length(n_vec))
  names(res_n) <- n_vec
  
  for (j in seq_along(n_vec)) {
    asy.n <- n_vec[j]
    
    ## check whether R0 itself gives a full-rank asymptotic covariance
    Psi_R0 <- metaSEM::asyCov(R0, n = asy.n)
    Psi_R0_ill <- (qr(Psi_R0)$rank != ncol(Psi_R0)) | (kappa(Psi_R0) > 1e8)
    
    if(!Psi_R0_ill){
      
      ## perturb around this R0
      Rhat <- rwaldcloud(R0, n = asy.n, B = B)
      
      Psi_Rhat <- lapply(1:B, function(i) {metaSEM::asyCov(x = Rhat[,,i], n = asy.n)})
      Psi_Rhat_ill <- sapply(Psi_Rhat, function(Q) {qr(Q)$rank != ncol(Q) | kappa(Q)> 1e8 })
      
      ## Wald distances: only compute for well-defined Psi_Rhat
      waldstat <- sapply(1:B, function(i) {
        if (Psi_Rhat_ill[i]) {
          NA_real_
        } else {
          Rsample <- Rhat[,,i]
          Rsample_vec <- Rsample[lower.tri(Rsample)]
          Psi_i <- Psi_Rhat[[i]]
          
          mahalanobis(
            x = Rsample_vec,
            center = R0_vec,
            cov = Psi_i
          )
        }
      })
      
      ## p-values
      pval <- pchisq(waldstat, df = df, lower.tail = FALSE)
      
      res_n[[j]] <- list(
        n = asy.n,
        waldstat = waldstat,
        pval = pval,
        ill = Psi_Rhat_ill,
        reverse_rate_all = mean(pval <= alpha),
        reverse_rate = mean(pval[!Psi_Rhat_ill]<= alpha, na.rm=TRUE),
        ill_sample_rate = mean(Psi_Rhat_ill)
      )
    } else {
      res_n[[j]] <- list(
        n = asy.n,
        waldstat = NA,
        pval = NA,
        ill = NA,
        reverse_rate_all = NA,
        reverse_rate = NA,
        ill_sample_rate = 1
      )
    }
  }
  
  res_all[[s]] <- list(
    R0 = R0,
    Psi_R0_ill = Psi_R0_ill,
    by_n = res_n
  )
}


#### Parallel Run
library(future.apply)
plan(multisession, workers = max(1, parallelly::availableCores() - 1))
run_one_s <- function(s, R0_all, n_vec, B, df, alpha) {
  cat("Running R0 system", s, "of", dim(R0_all)[1], "\n")
  
  ## generate one new R0
  R0 <- R0_all[s, , ]
  R0_vec <- R0[lower.tri(R0)]
  
  ## store results across n for this R0
  res_n <- vector("list", length(n_vec))
  names(res_n) <- n_vec
  
  for (j in seq_along(n_vec)) {
    asy.n <- n_vec[j]
    
    ## check whether R0 itself gives a full-rank asymptotic covariance
    Psi_R0 <- metaSEM::asyCov(R0, n = asy.n)
    Psi_R0_ill <- (qr(Psi_R0)$rank != ncol(Psi_R0)) | (kappa(Psi_R0) > 1e8)
    
    if (!Psi_R0_ill) {
      
      ## perturb around this R0
      Rhat <- rwaldcloud(R0, n = asy.n, B = B)
      
      Psi_Rhat <- lapply(1:B, function(i) {
        metaSEM::asyCov(x = Rhat[, , i], n = asy.n)
      })
      
      Psi_Rhat_ill <- sapply(Psi_Rhat, function(Q) {
        (qr(Q)$rank != ncol(Q)) | (kappa(Q) > 1e8)
      })
      
      ## Wald distances: only compute for well-defined Psi_Rhat
      waldstat <- sapply(1:B, function(i) {
        if (Psi_Rhat_ill[i]) {
          NA_real_
        } else {
          Rsample <- Rhat[, , i]
          Rsample_vec <- Rsample[lower.tri(Rsample)]
          Psi_i <- Psi_Rhat[[i]]
          
          mahalanobis(
            x = Rsample_vec,
            center = R0_vec,
            cov = Psi_i
          )
        }
      })
      
      ## p-values
      pval <- pchisq(waldstat, df = df, lower.tail = FALSE)
      
      res_n[[j]] <- list(
        n = asy.n,
        waldstat = waldstat,
        pval = pval,
        ill = Psi_Rhat_ill,
        reverse_rate_all = mean(pval <= alpha, na.rm = TRUE),
        reverse_rate = mean(pval[!Psi_Rhat_ill] <= alpha, na.rm = TRUE),
        ill_sample_rate = mean(Psi_Rhat_ill)
      )
      
    } else {
      res_n[[j]] <- list(
        n = asy.n,
        waldstat = rep(NA_real_, B),
        pval = rep(NA_real_, B),
        ill = rep(NA, B),
        reverse_rate_all = NA_real_,
        reverse_rate = NA_real_,
        ill_sample_rate = 1
      )
    }
  }
  
  list(
    R0 = R0,
    R0_well_defined = !all(sapply(res_n, function(x) isTRUE(x$ill_sample_rate == 1))),
    by_n = res_n
  )
}
res_all <- future_lapply(
  X = 1:S,
  FUN = run_one_s,
  R0_all = R0_all,
  n_vec = n_vec,
  B = B,
  df = df,
  alpha = alpha,
  future.seed = TRUE
)


## Summarize reverse containment rate
rate_mat <- sapply(1:S, function(s) {
  sapply(1:length(n_vec), function(j) {
    1-res_all[[s]]$by_n[[j]]$reverse_rate
  })
})
ill_mat <- sapply(1:S, function(s) {
  sapply(1:length(n_vec), function(j) {
    res_all[[s]]$by_n[[j]]$ill_sample_rate
  })
})

rownames(rate_mat) <- paste0("n=", n_vec)
rownames(ill_mat)  <- paste0("n=", n_vec)

rate_mat
ill_mat

## mean across R0 systems
rowMeans(rate_mat, na.rm=T)  # reverse containment rate
rowMeans(ill_mat, na.rm=T)

{
  kappa.i<- sapply(1:S, function(s) {
    P<- metaSEM::asyCov(R0_all[s,,], n = 1)
    kappa(P)
  })
  cols<- c("black", "blue","green")
  par(pin = c(3, 3), mar= c(9,4,4,2)+0.1)
  plot(log(kappa.i), rate_mat[1,], ylim=c(0,1.2), col = cols[1],pch = 16, cex = 1,
       xlab="log Kappa of Cov(R0)",
       ylab="Reverse Containment Rate", main =paste0("Reverse containment rate versus Condition Number of Cov(R0); p=",p))
  abline(h= 1-alpha, col="red", lty=2)
  points(log(kappa.i), rate_mat[2,], col = cols[2], pch = 16, cex = 1)
  points(log(kappa.i), rate_mat[3,], col = cols[3], pch = 16, cex = 1)
  legend("topleft", legend = paste("n =", n_vec, "; Average Rate = ",round(rowMeans(rate_mat, na.rm=T),2)), col = cols, pch = 16, bty = "n")
  mtext(paste0("R0 is sampled uniformly from the space of valid ",p," x ",p," correlation matrices. \n1000 Rhats are sampled from the Wald cloud centered at R0.\nReverse containment is assessed using the Mahalanobis distance between each Rhat and R0, using asyCov(Rhat).\nAverage reverse containment rate (excluding ill-conditioned asyCov(Rhat)) is reported for each R0"), 
        side=1, line=7, adj=0)
}

summary(log(kappa.i))

matplot(
  1:ncol(rate_mat),
  t(rate_mat),
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "R0 index",
  ylab = "Reverse Containment rate",
  ylim = c(0,1),
  main = "Reverse Containment rate for each R0"
)

legend(
  "topright",
  legend = paste0("n = ", n_vec),
  col = 1:length(n_vec),
  lty = 1,
  lwd = 2
)

abline(h = 1-alpha, col = "red", lty = 2)



##########################





out_dir = "/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/docs/plots/ReverseContainment2"
filename = paste0(c("nodes","ad","eta"), c(p, ad, eta), collapse="")
filename <- file.path(
  out_dir,
  paste0(filename, ".png")
)
#png(filename, width = 600, height= 500, res=300)


