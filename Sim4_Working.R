# PerturbR Sim3
rm(list=ls()); gc()
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2)
use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
py_config() 

source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R")  # Wald cloud
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R")  # Base CDA Functions

# Source python files
source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")

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



# Simulation --------------------------------------------------------------

# Full Simulation ---------------------------------------------------------

nTarget = 10
p= 6L; ad= 3L; asy.n= 2000L; ITER=10L
ps= p*(p-1)/2
alpha2<- 0.01   # used in CI test

perf_adj.waldcloud<- perf_adj.bootstrap<- perf_adj.unifcloud_02 <- perf_adj.unifcloud_01 <- perf_adj.unifcloud_005 <- perf_adj.unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
perf_ort.waldcloud<- perf_ort.bootstrap<- perf_ort.unifcloud_02 <- perf_ort.unifcloud_01 <- perf_ort.unifcloud_005 <- perf_ort.unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
perf_adj.boss_waldcloud <-perf_adj.boss_bootstrap <- perf_adj.boss_unifcloud_02 <- perf_adj.boss_unifcloud_01 <- perf_adj.boss_unifcloud_005 <- perf_adj.boss_unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
perf_ort.boss_waldcloud <-perf_ort.boss_bootstrap <- perf_ort.boss_unifcloud_02 <- perf_ort.boss_unifcloud_01 <- perf_ort.boss_unifcloud_005 <- perf_ort.boss_unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
# For objects above: First row is accuracy, Second row is sensitivity

# SIMULATION Starts here
{
  for(nT in 1:nTarget){
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ";Iteration:", nT, "\n")
    # Truth Layer:
    Target = er_dag_py(p=p, ad=ad, n= asy.n, K=1L)  # Consists of K data sets
    G0 <- Target$G
    R0 <- Target$R
    X<- Target$X[1,,]   # The only observed data
    
    Rhat <- cor(X) # Correspond to Step 3
    Rhat.waldcloud <- rwaldcloud(Rhat, n = asy.n, B = ITER )                  # Step 3a
    Rhat.unifcloud_001<- runifcloud(Rhat, n = asy.n, B = ITER, alpha =  0.01) # Step 3a
    #Rhat.unifcloud_005<- runifcloud(Rhat, n = asy.n, B = ITER, alpha =  0.05) # Step 3a
    Rhat.unifcloud_01<- runifcloud(Rhat, n = asy.n, B = ITER, alpha =  0.1)   # Step 3a
    #Rhat.unifcloud_02<- runifcloud(Rhat, n = asy.n, B = ITER, alpha =  0.2)   # Step 3a
    
    Rhat.bootstrap <- boot_cor(X, asy.n= asy.n, B = ITER)                     # Step 3b

    # Run CDA
    # PC
    {
      fit0<- pcalg::pc(suffStat= list(C = R0, n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2 )
      Ghat.0 <- as(fit0, "matrix")
      fit.waldcloud<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.waldcloud[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      #fit.corX<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.corX[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.unifcloud_001<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_001[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      #fit.unifcloud_005<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_005[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.unifcloud_01<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_01[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      #fit.unifcloud_02<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_02[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.bootstrap<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.bootstrap[,,i], n = asy.n/2), indepTest = gaussCItest, p= p, alpha = alpha2))
      
      Ghat.waldcloud<- lapply(fit.waldcloud, function(Y) as(Y, "matrix"))
      #Ghat.corX<- lapply(fit.corX, function(Y) as(Y, "matrix"))
      Ghat.unifcloud_001<- lapply(fit.unifcloud_001, function(Y) as(Y, "matrix"))
      #Ghat.unifcloud_005<- lapply(fit.unifcloud_005, function(Y) as(Y, "matrix"))
      Ghat.unifcloud_01<- lapply(fit.unifcloud_01, function(Y) as(Y, "matrix"))
      #Ghat.unifcloud_02<- lapply(fit.unifcloud_02, function(Y) as(Y, "matrix"))
      Ghat.bootstrap<- lapply(fit.bootstrap, function(Y) as(Y, "matrix"))
      
      # Compare Ghats with G0
      perf_adj.waldcloud[nT,,]<- sapply(Ghat.waldcloud, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      #perf.corX<- sapply(Ghat.corX, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric)) 
      perf_adj.unifcloud_001[nT,,]<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      #perf_adj.unifcloud_005[nT,,]<- sapply(Ghat.unifcloud_005, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_adj.unifcloud_01[nT,,]<- sapply(Ghat.unifcloud_01, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      #perf_adj.unifcloud_02[nT,,]<- sapply(Ghat.unifcloud_02, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_adj.bootstrap[nT,,]<- sapply(Ghat.bootstrap, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_ort.unifcloud_001[nT,,]<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      
      perf_ort.waldcloud[nT,,]<- sapply(Ghat.waldcloud, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.unifcloud_001[nT,,]<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      #perf_ort.unifcloud_005[nT,,]<- sapply(Ghat.unifcloud_005, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.unifcloud_01[nT,,]<- sapply(Ghat.unifcloud_01, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      #perf_ort.unifcloud_02[nT,,]<- sapply(Ghat.unifcloud_02, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.bootstrap[nT,,]<- sapply(Ghat.bootstrap, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    }
    
    
    # Boss
    {
      # Run Boss on truth R0
      Ghat_boss_0 <- boss_py(R0, asy.n)
      
      # Run Boss across clouds
      Ghat_boss_waldcloud <- lapply(1:ITER, function(i) boss_py(Rhat.waldcloud[,,i], asy.n))
      Ghat_boss_unifcloud_001 <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_001[,,i], asy.n))
      #Ghat_boss_unifcloud_005 <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_005[,,i], asy.n))
      Ghat_boss_unifcloud_01  <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_01[,,i],  asy.n))
      #Ghat_boss_unifcloud_02  <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_02[,,i],  asy.n))
      Ghat_boss_bootstrap     <- lapply(1:ITER, function(i) boss_py(Rhat.bootstrap[,,i], as.integer(asy.n/2)))
      
      # CDA metrics — adjacency
      perf_adj.boss_waldcloud[nT,,] <- sapply(Ghat_boss_waldcloud, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      perf_adj.boss_unifcloud_001[nT,,] <- sapply(Ghat_boss_unifcloud_001, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      #perf_adj.boss_unifcloud_005[nT,,] <- sapply(Ghat_boss_unifcloud_005, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      perf_adj.boss_unifcloud_01[nT,,] <- sapply(Ghat_boss_unifcloud_01, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      #perf_adj.boss_unifcloud_02[nT,,] <- sapply(Ghat_boss_unifcloud_02, function(ghat)  eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      perf_adj.boss_bootstrap[nT,,] <- sapply(Ghat_boss_bootstrap, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      
      # CDA metrics — orientation
      perf_ort.boss_waldcloud[nT,,] <- sapply(Ghat_boss_waldcloud, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      perf_ort.boss_unifcloud_001[nT,,] <- sapply(Ghat_boss_unifcloud_001, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      #perf_ort.boss_unifcloud_005[nT,,] <- sapply(Ghat_boss_unifcloud_005, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      perf_ort.boss_unifcloud_01[nT,,] <- sapply(Ghat_boss_unifcloud_01, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      #perf_ort.boss_unifcloud_02[nT,,] <- sapply(Ghat_boss_unifcloud_02, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      perf_ort.boss_bootstrap[nT,,] <- sapply(Ghat_boss_bootstrap, function(ghat) eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    }
  }
}


# ENSEMBLE
{
  criteria <- 0.5
  cvt <- function(mat, criteria) ifelse(mat>criteria, 1, 0 )
  
  Ghat.bootstrap.mean<- (Reduce(`+`, Ghat.bootstrap)/ITER)   %>% cvt(., criteria)
  Ghat.waldcloud.mean<- (Reduce(`+`, Ghat.waldcloud)/ITER)   %>% cvt(., criteria)
  Ghat.unifcloud_001.mean<- (Reduce(`+`, Ghat.unifcloud_001)/ITER)   %>% cvt(., criteria)
  Ghat.unifcloud_01.mean<- (Reduce(`+`, Ghat.unifcloud_01)/ITER)   %>% cvt(., criteria)
  
  Ghat.boss.bootstrap.mean<- (Reduce(`+`, Ghat_boss_bootstrap)/ITER)   %>% cvt(., criteria)
  Ghat.boss.waldcloud.mean<- (Reduce(`+`, Ghat_boss_waldcloud)/ITER)   %>% cvt(., criteria)
  Ghat.boss.unifcloud_001.mean<- (Reduce(`+`, Ghat_boss_unifcloud_001)/ITER)   %>% cvt(., criteria)
  Ghat.boss.unifcloud_01.mean<- (Reduce(`+`, Ghat_boss_unifcloud_01)/ITER)   %>% cvt(., criteria)
}

# Compare 
eval_cda(true_adj = G0, Ghat.bootstrap.mean)
eval_cda(true_adj = G0, Ghat.waldcloud.mean)
eval_cda(true_adj = G0, Ghat.unifcloud_001.mean)
eval_cda(true_adj = G0, Ghat.unifcloud_01.mean)


# Generate Plot
{
  perf_adj.list<- list(perf_adj.bootstrap, perf_adj.unifcloud_001, perf_adj.unifcloud_005, perf_adj.unifcloud_01, perf_adj.unifcloud_02)
  perf_ort.list<- list(perf_ort.bootstrap, perf_ort.unifcloud_001, perf_ort.unifcloud_005, perf_ort.unifcloud_01, perf_ort.unifcloud_02)
  perf_adj.bosslist<- list(perf_adj.boss_bootstrap, perf_adj.boss_unifcloud_001, perf_adj.boss_unifcloud_005, perf_adj.boss_unifcloud_01, perf_adj.boss_unifcloud_02)
  perf_ort.bosslist<- list(perf_ort.boss_bootstrap, perf_ort.boss_unifcloud_001, perf_ort.boss_unifcloud_005, perf_ort.boss_unifcloud_01, perf_ort.boss_unifcloud_02)
  
  perf_adj.pc<- lapply(perf_adj.list, function(dt) apply(dt,c(1,2),mean)) %>% simplify2array()
  perf_ort.pc<- lapply(perf_ort.list, function(dt) apply(dt,c(1,2),mean))%>% simplify2array()
  perf_adj.boss<- lapply(perf_adj.bosslist, function(dt) apply(dt,c(1,2),mean))%>% simplify2array()
  perf_ort.boss<- lapply(perf_ort.bosslist, function(dt) apply(dt,c(1,2),mean))%>% simplify2array()
  
  #rownames(perf_ort.boss)<-rownames(perf_adj.boss)<-rownames(perf_ort.pc)<-rownames(perf_adj.pc)<- c("Precision", "Sensitivity")
  #colnames(perf_ort.boss)<-colnames(perf_adj.boss)<-colnames(perf_ort.pc)<-colnames(perf_adj.pc)<- c("bootstrap", "unicloud_0.01","unicloud_0.05","unicloud_0.02","unicloud_0.2")
  
  #perf_adj.pc[,1,] # mean precision of adjacency for each graph
  # Distribution of mean performance across nTarget independently generated graphs
  mynames<- c("bootstrap", "99% Cloud","95% Cloud","90% Cloud","80% Cloud")
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  boxplot(perf_adj.pc[,1,], names= mynames, main="Adjacency Precision")
  boxplot(perf_adj.pc[,2,], names= mynames, main="Adjacency Sensitivity")
  boxplot(perf_ort.pc[,1,], names= mynames, main="Orientation Precision")
  boxplot(perf_ort.pc[,2,], names= mynames, main="Orientation Sensitivity")
  mtext("PC: Overall Performance Comparison", outer = TRUE, side = 3, line = 1.5, cex = 1.2)
  
  boxplot(perf_adj.boss[,1,], names= mynames, main="Adjacency Precision")
  boxplot(perf_adj.boss[,2,], names= mynames, main="Adjacency Sensitivity")
  boxplot(perf_ort.boss[,1,], names= mynames, main="Orientation Precision")
  boxplot(perf_ort.boss[,2,], names= mynames, main="Orientation Sensitivity")
  mtext("BOSS: Overall Performance Comparison", outer = TRUE, side = 3, line = 1.5, cex = 1.2)
  
}

perf_adj.pc
perf_ort.pc
perf_adj.boss
perf_ort.boss










# Full Simulation Codex ---------------------------------------------------

# ---- Settings (kept as requested) ----
nTarget <- 200L
p <- 10L
ad <- 2L
asy.n <- 20000L
ITER <- 300L
alpha2 <- 0.01


{
  methods <- c("waldcloud", "unifcloud_001", "unifcloud_01", "bootstrap")
  n_obs_by_method <- c(
    waldcloud = asy.n,
    unifcloud_001 = asy.n,
    unifcloud_01 = asy.n,
    bootstrap = as.integer(asy.n / 2L)
  )
  
  # perf arrays: [nT, metric(accuracy/sensitivity), iter]
  alloc_perf_array <- function() {
    array(
      NA_real_,
      dim = c(nTarget, 2L, ITER),
      dimnames = list(
        nT = as.character(seq_len(nTarget)),
        metric = c("accuracy", "sensitivity"),
        iter = as.character(seq_len(ITER))
      )
    )
  }
  
  perf <- list(
    pc = list(
      adj = setNames(lapply(methods, function(...) alloc_perf_array()), methods),
      orient = setNames(lapply(methods, function(...) alloc_perf_array()), methods)
    ),
    boss = list(
      adj = setNames(lapply(methods, function(...) alloc_perf_array()), methods),
      orient = setNames(lapply(methods, function(...) alloc_perf_array()), methods)
    )
  )
  
  # Store fitted adjacency matrices as arrays [p, p, ITER] for each nT and method
  fit <- list(
    truth    = vector("list", nTarget),
    baseline = list(pc = vector("list", nTarget), boss = vector("list", nTarget)),
    optimal  = list(pc = vector("list", nTarget), boss = vector("list", nTarget))
  )
  
  # fit[[method]][[nT]]$algo[[i]]  (i = 1..ITER)
  for (m in methods) {
    fit[[m]] <- vector("list", nTarget)
    for (nT in seq_len(nTarget)) {
      fit[[m]][[nT]] <- list(
        pc   = vector("list", ITER),
        boss = vector("list", ITER)
      )
    }
  }
  
  
  run_pc_mat <- function(Cmat, n_obs) {
    as(
      pcalg::pc(
        suffStat = list(C = Cmat, n = n_obs),
        indepTest = gaussCItest,
        p = p,
        alpha = alpha2
      ),
      "matrix"
    )
  }
  
  get_metric <- function(G0, Ghat) {
    out <- eval_cda(true_adj = t(G0), est_adj = t(Ghat))
    list(adj = out$metric.adj, orient = out$metric.orient)
  }
  
}

for (nT in seq_len(nTarget)) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "; Iteration:", nT, "\n")
  
  Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 1L)
  G0 <- Target$G
  R0 <- Target$R
  
  X <- Target$X[1, , ]
  Rhat <- cor(X)
  
  Rhat_arr <- list(
    waldcloud     = rwaldcloud(Rhat, n = asy.n, B = ITER),
    unifcloud_001 = runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.01),
    unifcloud_01  = runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.1),
    bootstrap     = boot_cor(X, asy.n = asy.n, B = ITER)
  )
  
  fit$truth[[nT]] <- G0
  
  # If "optimal" truly means using R0 (true correlation), use R0 here
  fit$optimal$pc[[nT]]   <- run_pc_mat(R0, asy.n)
  fit$optimal$boss[[nT]] <- boss_py(R0, asy.n)
  
  fit$baseline$pc[[nT]]   <- run_pc_mat(Rhat, asy.n)
  fit$baseline$boss[[nT]] <- boss_py(Rhat, asy.n)
  
  for (m in methods) {
    n_obs <- n_obs_by_method[[m]]
    
    for (i in seq_len(ITER)) {
      fit[[m]][[nT]]$pc[[i]]   <- run_pc_mat(Rhat_arr[[m]][, , i], n_obs)
      fit[[m]][[nT]]$boss[[i]] <- boss_py(  Rhat_arr[[m]][, , i], n_obs)
    }
  }
}

#fit$truth: G0
#fit$optimal: learn based on R0
#fit$baseline: learn based on Rhat
#fit$waldcloud and other methods.... learn based on perturbation of Rhat

{ #Build Comparison Table
  library(data.table)
  # out: 4x3, rows TP/FP/FN/TN, col2 adj, col3 orient
  out_to_row <- function(out_mat) {
    data.table(
      TP_adj = as.integer(out_mat[1, 2]),
      FP_adj = as.integer(out_mat[2, 2]),
      FN_adj = as.integer(out_mat[3, 2]),
      TN_adj = as.integer(out_mat[4, 2]),
      TP_ort = as.integer(out_mat[1, 3]),
      FP_ort = as.integer(out_mat[2, 3]),
      FN_ort = as.integer(out_mat[3, 3]),
      TN_ort = as.integer(out_mat[4, 3])
    )
  }
  
  build_eval_table <- function(fit,
                               methods = c("baseline","waldcloud","unifcloud_001","unifcloud_01","bootstrap"),
                               algos   = c("pc","boss"),
                               targets = c("G0","G0hat","Ghat"),
                               nTarget,
                               ITER) {
    
    # allocate generously; baseline has only 1 i, others have ITER
    max_rows <- nTarget * length(algos) * length(targets) * (
      1L + (length(methods) - 1L) * ITER
    )
    
    rows <- vector("list", max_rows)
    idx <- 0L
    
    for (algo in algos) {
      for (nT in seq_len(nTarget)) {
        
        # Targets for this replicate (G0 shared; others algo-specific)
        G0    <- fit$truth[[nT]]
        G0hat <- fit$optimal[[algo]][[nT]]
        Ghat  <- fit$baseline[[algo]][[nT]]
        
        for (m in methods) {
          
          # baseline has a single estimate (i=1); perturb methods have i=1..ITER
          i_seq <- if (m == "baseline") 1L else seq_len(ITER)
          
          for (i in i_seq) {
            
            # estimator graph
            G_est <- if (m == "baseline") {
              fit$baseline[[algo]][[nT]]
            } else {
              fit[[m]][[nT]][[algo]][[i]]
            }
            
            for (tgt in targets) {
              
              G_true <- switch(
                tgt,
                G0    = G0,
                G0hat = G0hat,
                Ghat  = Ghat,
                stop("Unknown target: ", tgt)
              )
              
              out_mat <- eval_cda(true_adj = as.matrix(G_true),
                                  est_adj  = as.matrix(G_est))$out
              
              idx <- idx + 1L
              rows[[idx]] <- cbind(
                data.table(method = m, algo = algo, target = tgt, nT = nT, i = i),
                out_to_row(out_mat)
              )
            }
          }
        }
      }
    }
    
    rbindlist(rows[seq_len(idx)], use.names = TRUE)
  }
  
  dt_eval <- build_eval_table(
    fit = fit,
    methods = c("baseline", "waldcloud","unifcloud_001","unifcloud_01","bootstrap"),
    algos   = c("pc","boss"),
    targets = c("G0","G0hat","Ghat"),
    nTarget = nTarget,
    ITER    = ITER
  )
  
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  dt_eval[, `:=`(
    precision_adj = safe_div(TP_adj, TP_adj + FP_adj),
    recall_adj    = safe_div(TP_adj, TP_adj + FN_adj),
    precision_ort = safe_div(TP_ort, TP_ort + FP_ort),
    recall_ort    = safe_div(TP_ort, TP_ort + FN_ort)
  )][, `:=`(
    f1_adj        = safe_div(2 * precision_adj * recall_adj, precision_adj + recall_adj),
    f1_ort        = safe_div(2 * precision_ort * recall_ort, precision_ort + recall_ort)
  )]

}


dt_eval[method=="baseline"]

{ # Plot
  dt_plot <- melt(
    dt_eval,
    id.vars = c("method","algo","target","nT","i"),
    measure.vars = c("precision_adj","recall_adj","f1_adj",
                     "precision_ort","recall_ort","f1_ort"),
    variable.name = "measure",
    value.name = "score"
  )
  
  dt_plot[, eval_type := ifelse(grepl("_adj$", measure), "adj", "orient")]
  dt_plot[, metric := sub("_(adj|ort)$", "", measure)]
  dt_plot[, measure := NULL]
  dt_nt_median <- dt_plot[
    , .(score_med = median(score, na.rm = TRUE)),
    by = .(method, algo, target, nT, eval_type, metric)
  ]
  
  plot_nt_median_box <- function(dt_nt_median,
                                 algo_pick = "pc",
                                 eval_type_pick = "adj",
                                 metric_pick = "precision",
                                 method_levels = c("baseline","waldcloud","unifcloud_001","unifcloud_01","bootstrap"),
                                 target_levels = c("G0","G0hat","Ghat")) {
    
    df <- dt_nt_median[
      algo == algo_pick &
        eval_type == eval_type_pick &
        metric == metric_pick
    ]
    
    df[, method := factor(method, levels = method_levels)]
    df[, target := factor(target, levels = target_levels)]
    
    ggplot(df, aes(x = method, y = score_med, fill = target)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.7) +
      labs(
        x = "Method",
        y = paste0("Median ", metric_pick, " (per nT)"),
        title = paste0(toupper(algo_pick), " — ", eval_type_pick, " ", metric_pick, " (median over i)"),
        fill = "Compared to"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }
  
  
}
par(mfrow=c(2,3))
plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="adj", metric_pick="precision")
plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="adj", metric_pick="recall")
plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="adj", metric_pick="f1")
plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="orient", metric_pick="precision")
plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="orient", metric_pick="recall")
plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="orient", metric_pick="f1")

plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="adj", metric_pick="precision")
plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="adj", metric_pick="recall")
plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="adj", metric_pick="f1")
plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="orient", metric_pick="precision")
plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="orient", metric_pick="recall")
plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="orient", metric_pick="f1")





# Ensemble results
k <- 0.6  # proportion threshold
methods <- c("waldcloud", "unifcloud_001", "unifcloud_01", "bootstrap")
ensemble     <- setNames(vector("list", length(methods)), methods)
ensemble.bin <- setNames(vector("list", length(methods)), methods)

for (m in methods) {
  ensemble[[m]]     <- vector("list", nTarget)
  ensemble.bin[[m]] <- vector("list", nTarget)
  
  for (nT in seq_len(nTarget)) {
    
    pc_sum   <- Reduce(`+`, fit[[m]][[nT]]$pc[seq_len(ITER)])
    boss_sum <- Reduce(`+`, fit[[m]][[nT]]$boss[seq_len(ITER)])
    
    # store counts (0..ITER)
    ensemble[[m]][[nT]] <- list(pc = pc_sum, boss = boss_sum)
    
    # threshold on frequency
    ensemble.bin[[m]][[nT]] <- list(
      pc   = ifelse((pc_sum   / ITER) >= k, 1L, 0L),
      boss = ifelse((boss_sum / ITER) >= k, 1L, 0L)
    )
  }
}

ensemble.bin$waldcloud[[1]]$pc

algorithms <- names(ensemble.bin)
compare_res <- setNames(vector("list", length(algorithms)), algorithms)
for (algo in algorithms) {
  methods <- names(ensemble.bin[[algo]][[1]])
  compare_res[[algo]] <- setNames(vector("list", length(methods)), methods)
  
  for (m in methods) {
    ensemble.adj <- matrix(NA_real_, nrow = nTarget, ncol = 4L)
    ensemble.ort <- matrix(NA_real_, nrow = nTarget, ncol = 4L)
    
    for (i in seq_len(nTarget)) {
      truth_i <- fit$truth[[i]]  
      truth_i <- fit$point[[algo]][[i]]   
      est_i   <- ensemble.bin[[algo]][[i]][[m]]
      
      colnames(truth_i)<- rownames(truth_i)<- NULL
      colnames(est_i)<- rownames(est_i)<- NULL
      
      class(truth_i)
      class(est_i)
      
      out <- eval_cda(true_adj = as.matrix(truth_i), est_adj = as.matrix(est_i))$out
      ensemble.adj[i, ] <- out[, 2]
      ensemble.ort[i, ] <- out[, 3]
    }
    
    colnames(ensemble.adj) <- c("TP", "FP", "FN", "TN")
    colnames(ensemble.ort) <- c("TP", "FP", "FN", "TN")
    
    compare_res[[algo]][[m]] <- list(
      adj = ensemble.adj,
      ort = ensemble.ort
    )
  }
}


safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
add_prf <- function(conf_mat) {
  # conf_mat has columns: TP, FP, FN, TN
  tp <- conf_mat[, "TP"]
  fp <- conf_mat[, "FP"]
  fn <- conf_mat[, "FN"]
  
  precision <- safe_div(tp, tp + fp)
  recall <- safe_div(tp, tp + fn)
  f1 <- safe_div(2 * precision * recall, precision + recall)
  
  cbind(conf_mat, precision = precision, recall = recall, f1 = f1)
}

# Compute metrics for both adjacency and orientation
compare_metrics <- lapply(compare_res, function(by_method) {
  lapply(by_method, function(x) {
    list(
      adj = add_prf(x$adj),
      ort = add_prf(x$ort)
    )
  })
})
compare_summary <- lapply(compare_metrics, function(by_method) { # summary by each method across all nT
  lapply(by_method, function(x) {
    list(
      adj = apply(x$adj[, c("precision", "recall", "f1")],2, quantile, c( 0.5 ), na.rm = TRUE),
      ort = apply(x$ort[, c("precision", "recall", "f1")],2, quantile, c( 0.5), na.rm = TRUE)
    )
  })
})

summary_tables <- lapply(names(compare_summary), function(algo) {
  by_method <- compare_summary[[algo]]
  
  adj_df <- do.call(
    rbind,
    lapply(names(by_method), function(method) {
      x <- by_method[[method]]$adj[c("precision", "recall", "f1")]
      data.frame(
        method = method,
        precision = unname(x["precision"]),
        recall = unname(x["recall"]),
        f1 = unname(x["f1"]),
        row.names = NULL
      )
    })
  )
  
  ort_df <- do.call(
    rbind,
    lapply(names(by_method), function(method) {
      x <- by_method[[method]]$ort[c("precision", "recall", "f1")]
      data.frame(
        method = method,
        precision = unname(x["precision"]),
        recall = unname(x["recall"]),
        f1 = unname(x["f1"]),
        row.names = NULL
      )
    })
  )
  
  rownames(adj_df) <- adj_df$method
  adj_df$method <- NULL
  
  rownames(ort_df) <- ort_df$method
  ort_df$method <- NULL
  
  list(adj = adj_df, ort = ort_df)
})
names(summary_tables) <- names(compare_summary)
summary_tables



# CALIBRATION







fit$truth
fit_sum$pc[[1]]$waldcloud 


# Microbenchmark ----------------------------------------------------------
library(microbenchmark)
p= 15L; ad= 3L; asy.n= 20000L; ITER=300L
ps= p*(p-1)/2
#alpha1<- 0.05   # alpha1= 1 minus coverage level
alpha2<- 0.01   # used in CI test

{
  # Truth Layer:
  Target = er_dag_py(p=p, ad=ad, n= asy.n, K=1L)  # Consists of K data sets
  G0 <- Target$G
  R0 <- Target$R
  BETA <- Target$B
  OMEGA <- Target$O
  X<- Target$X[1,,]   # The only observed data
  Rhat<- cor(X)
}


{
  #MICROBENCHMARK
  res<- microbenchmark::microbenchmark(
    "BOOTSTRAP"= {boot_cor(X= X, asy.n= asy.n, B= ITER)},
    "FASTR" = {runifcloud(R0=Rhat, n = asy.n, B = ITER, alpha =  0.01)},
    times=10
  )
}
res


library(microbenchmark)

# parameter grid
ad   <- 3L
ITER <- 300L
times<- 100

grid <- expand.grid(
  p      = c(5L, 10L, 20L, 50L),
  asy.n  = c(2000L, 5000L, 20000L, 50000L),
  alpha2 = c(0.10, 0.01),
  KEEP.OUT.ATTRS = FALSE
)

bench_one <- function(p, ad, asy.n, ITER, alpha2, times = 10) {
  # --- setup OUTSIDE benchmark ---
  Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 1L)
  X      <- Target$X[1,,]
  Rhat   <- cor(X)
  
  # --- benchmark ONLY the calls ---
  res <- microbenchmark(
    BOOTSTRAP = boot_cor(X = X, asy.n = asy.n, B = ITER),
    FASTR     = runifcloud(R0 = Rhat, n = asy.n, B = ITER, alpha = alpha2),
    times = times
  )
  
  # return a tidy data.frame with params attached
  data.frame(
    p = p, asy.n = asy.n, alpha2 = alpha2,
    expr = res$expr,
    time_ns = res$time,
    stringsAsFactors = FALSE
  )
}



all_res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  with(grid[i, ], bench_one(p = p, ad = ad, asy.n = asy.n, ITER = ITER, alpha2 = alpha2, times= times))
}))

# summarize (median ms) by config + method
summary_res <- aggregate(
  time_ns ~ p + asy.n + alpha2 + expr,
  data = all_res,
  FUN = function(x) median(x) / 1e6
)
names(summary_res)[names(summary_res) == "time_ns"] <- "median_ms"

summary_res[order(summary_res$p, summary_res$asy.n, summary_res$alpha2, summary_res$expr), ]





# bench::mark -------------------------------------------------------------
pacman::p_load(bench)
bench_one <- function(p, ad, asy.n, ITER, alpha2,
                      times = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- setup OUTSIDE benchmarking ---
  Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 1L)
  X      <- Target$X[1,,]
  Rhat   <- cor(X)
  
  # --- benchmark ONLY the calls ---
  bm <- bench::mark(
    BOOTSTRAP = boot_cor(X = X, asy.n = asy.n, B = ITER),
    FASTR     = runifcloud(R0 = Rhat, n = asy.n, B = ITER, alpha = alpha2),
    iterations = times,
    check = FALSE
  )
  
  # --- tidy output (one row per method) ---
  
  out<- data.table(
    p = p,
    asy.n = asy.n,
    alpha2 = alpha2,
    as.data.table(bm)
  )
  out[,expression:= names(bm$expression)]
  
  return(out)
  
}

# parameter grid
ad   <- 3L
ITER <- 300L
times<- 100

grid <- expand.grid(
  p      = c(5L, 10L, 15L, 20L),
  asy.n  = c(2000L, 5000L, 20000L, 50000L),
  alpha2 = c(0.10, 0.01),
  KEEP.OUT.ATTRS = FALSE
)

all_bench <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  with(grid[i, ], bench_one(
    p = p, ad = ad, asy.n = asy.n, ITER = ITER, alpha2 = alpha2,
    times = times,
    seed = 1000 + i          # optional: stabilize across grid points
  ))
}))

all_bench[,`:=`(memory=NULL, time=NULL, gc=NULL )]
colnames(all_bench)


library(ggplot2)
ggplot(all_bench,
       aes(x = p, y = median, color = expression, group = expression)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(alpha2 ~ asy.n, labeller = label_both) +
  labs(
    x = "Number of variables (p)",
    y = "Median runtime (ms)",
    color = "Method",
    title = "Runtime Comparison"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
ggplot(all_bench,
       aes(x = p, y = mem_alloc, color = expression, group = expression)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(alpha2 ~ asy.n, labeller = label_both) +
  labs(
    x = "Number of variables (p)",
    y = "Memory allocation",
    color = "Method",
    title = "Memory allocation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
ggplot(all_bench,
       aes(x = p, y = n_gc, color = expression, group = expression)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(alpha2 ~ asy.n, labeller = label_both) +
  labs(
    x = "Number of variables (p)",
    y = "n_gc",
    color = "Method",
    title = "Total number of garbage collections performed over all iterations"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
#saveRDS(all_bench,"~/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/plots/CompareComputingResources/all_bench.Rds" )

BOOTSTRAP = boot_cor(X = X, asy.n = asy.n, B = ITER)
FASTR     = runifcloud(R0 = Rhat, n = asy.n, B = ITER, alpha = alpha2)
object.size(BOOTSTRAP)
object.size(FASTR)
60224/(1024)
all_bench$mem_alloc %>% str()
all_bench[1:2,]


Rprofmem("mem.out")
test1<- boot_cor(X = X, asy.n = asy.n, B = ITER)
Rprofmem(NULL)
mem_lines<- readLines("mem.out")
alloc_bytes <- as.numeric(sub(" .*", "", mem_lines))
alloc_bytes <- alloc_bytes[!is.na(alloc_bytes)]
sum(alloc_bytes)/(1024^2)


Rprofmem("mem.out")
test2<- runifcloud(R0 = Rhat, n = asy.n, B = ITER, alpha = alpha2)
Rprofmem(NULL)
mem_lines<- readLines("mem.out")
alloc_bytes <- as.numeric(sub(" .*", "", mem_lines))
alloc_bytes <- alloc_bytes[!is.na(alloc_bytes)]
sum(alloc_bytes)/(1024^2)

15*15*300*8/1024^2
