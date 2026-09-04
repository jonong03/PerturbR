# Import functions --------------------------------------------------------
pacman::p_load(xtable, future.apply, data.table, rethinking, dplyr, reticulate, mvtnorm, Rcpp)
source("https://raw.githubusercontent.com/jonong03/PerturbR/main/docs/WorkingCode/Functions.R")
cpp_url <- "https://raw.githubusercontent.com/jonong03/PerturbR/main/docs/WorkingCode/calc_psi.cpp"
tmp <- tempfile(fileext = ".cpp")
download.file(cpp_url, tmp)
Rcpp::sourceCpp(tmp)

#use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
#py_config() 
#source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
#source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")

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



# Parameter Matrix --------------------------------------------------------
{
  gc()
  R2true <- c(0.5)
  phi <- c(0, 0.2, 0.5, 0.8)
  ss <- c(1000, 3000, 5000)
  p <- c(2)
  
  beta.list <- c(
    "function(p) rep(0.6, p)",
    "function(p) rep(0.2, p)",
    "function(p) rep(c(-0.4, 0.4), length= p)"
  )
  
  param.m <- expand.grid(R2true = R2true, phi = phi, ss = ss, p = p, beta = beta.list,stringsAsFactors = FALSE) %>% as.data.table
  param.m$beta_value <- Map(
    function(f, p) {
      fun <- eval(parse(text = f))
      fun(p)
    },
    param.m$beta,
    param.m$p
  )
  #param.m[, scenarioID := .I]
  #param.m <- param.m[rep(seq_len(nrow(param.m)), each = 3), ]
  #param.m[, method := rep(c("analytical", "bootstrap", "perturbR"), length.out = .N)]
}
{
  ncores <- detectCores()-1
  B = 500 # Number of Monte Carlo Replicates
  BOOT.ITER = 300   #J
  nchain = 3000
}
param.m


# Data Generation ---------------------------------------------------------
{
  param.m[, DT.list := Map(
    datagen,
    p, beta_value, ss, phi, R2true, ITER= B)]
}

# Analytical CI -----------------------------------------------------------

FL_eval <- function(DT.list, p, ss, R2true) {
  
  # Rhat = cor(DT), then compute R2
  R2hat <- sapply(DT.list, function(DT) computeRsq(cor(as.matrix(DT))))
  
  # CI for each Monte Carlo replicate
  fl.ci <- sapply(R2hat, function(r2) MBESS::ci.R2(R2 = r2, K = p, N = ss))
  lb <- unlist(fl.ci[1, ])
  ub <- unlist(fl.ci[3, ])
  
  c(coverage.fl = mean(R2true >= lb & R2true <= ub),
    width.fl    = mean(ub - lb))
}
fl.result <- Map(
  FL_eval,
  param.m$DT.list, param.m$p, param.m$ss, param.m$R2true)


# Bootstrap CI ------------------------------------------------------------
BOOT_eval <- function(DT.list, R2true, BOOT.ITER) {
  boot.ci <- lapply(DT.list, function(DT) {
    empR2dist <- sapply(seq_len(BOOT.ITER), function(i) {
      boot.id <- sample.int(nrow(DT), nrow(DT), replace = TRUE)
      DTboot <- DT[boot.id, , drop = FALSE]
      computeRsq(cor(DTboot))
    })
    quantile(empR2dist, c(0.025, 0.975))
  })
  
  boot.ci <- do.call(cbind, boot.ci)
  
  c(
    coverage.boot = mean(R2true >= boot.ci[1, ] & R2true <= boot.ci[2, ]),
    width.boot = mean(boot.ci[2, ] - boot.ci[1, ])
  )
}
boot.result <- parallel::mclapply(
  seq_len(nrow(param.m)),
  function(i) {
    tryCatch(
      BOOT_eval(DT.list = param.m$DT.list[[i]],R2true = param.m$R2true[i],BOOT.ITER = BOOT.ITER),
      error = function(e) {
        list(
          scenarioID = param.m$scenarioID[i],
          row = i,
          error = conditionMessage(e)
        )
      }
    )
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)
boot.result <- do.call(rbind, boot.result)
boot.result

# PerturbR ----------------------------------------------------------------
perturbR_eval <- function(DT.list, ss, R2true, nchain) {
  
  perturbR.ci <- lapply(DT.list, function(DT) {
    theta <- cor(DT)
    thetacloud <- make_mcmc_fast(S = theta, asy.n = ss, NCHAIN = nchain, alpha = 0.05)
    empR2dist <- vapply(seq_len(nchain), function(i) computeRsq(thetacloud[,,i]),numeric(1))
    quantile(empR2dist, c(0, 1))
  })
  
  perturbR.ci <- do.call(cbind, perturbR.ci)
  
  c(
    coverage.perturbR = mean(R2true >= perturbR.ci[1, ] & R2true <= perturbR.ci[2, ]),
    width.perturbR = mean(perturbR.ci[2, ] - perturbR.ci[1, ])
  )
}
perturbR.result <- parallel::mclapply(
  seq_len(nrow(param.m)),
  function(i) {
    
    tryCatch(
      perturbR_eval(
        DT.list = param.m$DT.list[[i]],
        ss = param.m$ss[i],
        R2true = param.m$R2true[i],
        nchain = nchain
      ),
      error = function(e) {
        list(
          scenarioID = param.m$scenarioID[i],
          row = i,
          error = conditionMessage(e)
        )
      }
    )
    
  },
  mc.cores = ncores,
  mc.preschedule = FALSE
)
perturbR.result <- do.call(rbind, perturbR.result)



# Assemble Results --------------------------------------------------------

out1<- copy(param.m[,.(scenarioID, R2true, phi, ss, p, beta_value)])
out1[,`:=`(
  coverage.fl = as.vector(sapply(fl.result, `[`, "coverage.fl")),
  width.fl    = as.vector(sapply(fl.result, `[`, "width.fl")),
  coverage.boot = boot.result[, "coverage.boot"],
  width.boot    = boot.result[, "width.boot"],
  coverage.perturbR = perturbR.result[, "coverage.perturbR"],
  width.perturbR    = perturbR.result[, "width.perturbR"]
)]
out1

fout <- paste0("UseCase/msi/primary_out1.rds")
saveRDS(out1, file = fout)