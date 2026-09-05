# ========================================================================
# PerturbR R2 Simulation - Primary Analysis
# One SLURM array task = one parameter scenario
# ========================================================================


# Packages ----------------------------------------------------------------

library(data.table)
library(mvtnorm)
library(Rcpp)
library(MBESS)

source("docs/WorkingCode/Functions.R")
Rcpp::sourceCpp("docs/WorkingCode/calc_psi.cpp")
source("UseCase/msi/params.R")


# Read SLURM task ID -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("No scenario ID supplied.")
}

task_id <- as.integer(args[1])

if (is.na(task_id)) {
  stop("Scenario ID must be an integer.")
}

if (task_id < 1 || task_id > nrow(param.m)) {
  stop(sprintf("Invalid scenario ID %d. param.m contains %d scenarios.",
               task_id, nrow(param.m)))
}

param.m <- param.m[scenarioID == task_id]

if (nrow(param.m) != 1) {
  stop("Expected exactly one parameter scenario.")
}

cat("Running scenario:", task_id, "\n")
print(param.m)


# Simulation settings -----------------------------------------------------

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

B <- 500
BOOT.ITER <- 300
nchain <- 3000

cat("Monte Carlo replicates:", B, "\n")
cat("Bootstrap samples:", BOOT.ITER, "\n")
cat("PerturbR samples:", nchain, "\n")
cat("CPUs:", ncores, "\n")


# Reproducibility ---------------------------------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(20260406 + task_id)


# Helper functions --------------------------------------------------------

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
  datagen <- function(p, beta_value, ss, phi, R2true, ITER) {
    beta <- matrix(beta_value, nrow = p, ncol = 1)
    Rxx <- gen_AR1(p = p, rho = phi)
    S <- t(beta) %*% Rxx %*% beta
    sigma2 <- as.numeric(S) * (1 - R2true) / R2true
    
    n <- ss * ITER
    X <- rmvnorm(n=n, mean= rep(0,p), sigma= Rxx)
    eps <- rnorm(n=n, mean= 0, sd = sqrt(sigma2))
    
    Y <- X %*% beta + eps
    
    # Indices corresponding to each Monte Carlo replicate
    idx <- split(seq_len(n), rep(seq_len(ITER), each = ss))
    
    # Split into datasets
    DT.list <- lapply(
      idx,
      function(ii) as.matrix(cbind(Y = Y[ii], X[ii, , drop = FALSE]))
    )
    
    return(DT.list)
  }
  computeRsq_beta <- function(beta, theta) {
    2 * sum(beta * theta[-1, 1]) -
      drop(beta %*% theta[-1, -1] %*% beta)
  }
  split_DT <- function(DT.list, train_prop = 0.7) {
    
    lapply(DT.list, function(DT) {
      
      n <- nrow(DT)
      train_id <- sample.int(n, size = floor(train_prop * n))
      
      list(
        train = DT[train_id, , drop = FALSE],
        test  = DT[-train_id, , drop = FALSE]
      )
    })
  }
  estimate_beta <- function(theta) {
    theta_xx <- theta[-1, -1]
    theta_xy <- theta[-1, 1]
    beta_hat <- solve(theta_xx, theta_xy)
    as.numeric(beta_hat)
  }
}

# Data generation ---------------------------------------------------------

cat("\nGenerating data...\n")

DT.list <- datagen(
  p = param.m$p[1],
  beta_value = param.m$beta_value[[1]],
  ss = param.m$ss[1],
  phi = param.m$phi[1],
  R2true = param.m$R2true[1],
  ITER = B
)

cat("Generated", length(DT.list), "Monte Carlo datasets.\n")


# Analytical CI -----------------------------------------------------------

FL_eval <- function(DT.list, p, ss, R2true) {
  
  R2hat <- vapply(
    DT.list,
    function(DT) computeRsq(cor(DT)),
    numeric(1)
  )
  
  fl.ci <- lapply(
    R2hat,
    function(r2) MBESS::ci.R2(R2 = r2, K = p, N = ss)
  )
  
  lb <- vapply(
    fl.ci,
    function(x) x[["Lower.Conf.Limit.R2"]],
    numeric(1)
  )
  
  ub <- vapply(
    fl.ci,
    function(x) x[["Upper.Conf.Limit.R2"]],
    numeric(1)
  )
  
  return(c(
    coverage.fl = mean(R2true >= lb & R2true <= ub),
    width.fl = mean(ub - lb)
  ))
}

cat("\nRunning analytical CI...\n")

fl.result <- FL_eval(
  DT.list = DT.list,
  p = param.m$p[1],
  ss = param.m$ss[1],
  R2true = param.m$R2true[1]
)

print(fl.result)


# Bootstrap CI ------------------------------------------------------------

BOOT_eval <- function(DT.list, R2true, BOOT.ITER, ncores) {
  
  cores <- min(ncores, length(DT.list))
  
  boot.ci <- parallel::mclapply(
    DT.list,
    function(DT) {
      
      empR2dist <- vapply(
        seq_len(BOOT.ITER),
        function(i) {
          boot.id <- sample.int(nrow(DT), size = nrow(DT), replace = TRUE)
          DTboot <- DT[boot.id, , drop = FALSE]
          
          computeRsq(cor(DTboot))
        },
        numeric(1)
      )
      
      quantile(empR2dist, probs = c(0.025, 0.975), names = FALSE)
    },
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = TRUE
  )
  
  boot.ci <- do.call(cbind, boot.ci)
  
  return(c(
    coverage.boot = mean(R2true >= boot.ci[1, ] & R2true <= boot.ci[2, ]),
    width.boot = mean(boot.ci[2, ] - boot.ci[1, ])
  ))
}


cat("\nRunning bootstrap CI...\n")

boot.result <- BOOT_eval(
  DT.list = DT.list,
  R2true = param.m$R2true[1],
  BOOT.ITER = BOOT.ITER,
  ncores = ncores
)

print(boot.result)


# PerturbR CI -------------------------------------------------------------

perturbR_eval <- function(DT.list, ss, R2true, nchain, ncores) {
  
  cores <- min(ncores, length(DT.list))
  
  perturbR.ci <- parallel::mclapply(
    DT.list,
    function(DT) {
      
      theta <- cor(DT)
      
      thetacloud <- make_mcmc_fast(
        S = theta,
        asy.n = ss,
        NCHAIN = nchain,
        alpha = 0.05
      )
      
      empR2dist <- vapply(
        seq_len(nchain),
        function(i) computeRsq(thetacloud[, , i]),
        numeric(1)
      )
      
      range(empR2dist)
    },
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = TRUE
  )
  
  perturbR.ci <- do.call(cbind, perturbR.ci)
  
  return(c(
    coverage.perturbR = mean(R2true >= perturbR.ci[1, ] &
                               R2true <= perturbR.ci[2, ]),
    width.perturbR = mean(perturbR.ci[2, ] - perturbR.ci[1, ])
  ))
}


cat("\nRunning PerturbR CI...\n")

perturbR.result <- perturbR_eval(
  DT.list = DT.list,
  ss = param.m$ss[1],
  R2true = param.m$R2true[1],
  nchain = nchain,
  ncores = ncores
)

print(perturbR.result)


# Assemble results --------------------------------------------------------

out1 <- copy(
  param.m[, .(scenarioID, R2true, phi, ss, p, beta, beta_value)]
)

out1[, `:=`(
  coverage.fl = unname(fl.result["coverage.fl"]),
  width.fl = unname(fl.result["width.fl"]),
  coverage.boot = unname(boot.result["coverage.boot"]),
  width.boot = unname(boot.result["width.boot"]),
  coverage.perturbR = unname(perturbR.result["coverage.perturbR"]),
  width.perturbR = unname(perturbR.result["width.perturbR"])
)]

print(out1)


# Save results ------------------------------------------------------------

job_id <- Sys.getenv("SLURM_ARRAY_JOB_ID")

job_dir <- file.path(
  "UseCase/msi/results",
  paste0("job_", job_id)
)

dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

fout <- file.path(
  job_dir,
  sprintf("primary_scenario_%04d.rds", task_id)
)

saveRDS(out1, file = fout)

cat("\nJob ID:", job_id, "\n")
cat("Scenario:", task_id, "\n")
cat("Saved:", fout, "\n")