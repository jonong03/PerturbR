# Use case 1: CDA

{
  rm(list=ls()); gc()
  pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr, plotly)
  plan(sequential)
  
  use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
  py_config() 
  reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
  reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
  reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py_cpdag.py")
  
  source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
  source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)
  make_boot <- function(Xi, ITER) {
    ss_i <- nrow(Xi)
    p <- ncol(Xi)
    out <- array(NA_real_, dim = c(p, p, ITER))
    
    for (b in seq_len(ITER)) {
      idx <- sample(ss_i, ss_i, replace = TRUE)
      out[, , b] <- cor(Xi[idx,])
    }
    out
  }
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
      Psi_Rcurrent<- metaSEM::asyCov(Rcurrent, n = asy.n)
      
      Rs<- rwaldcloud(Rcurrent, Psi= scale_factor*Psi_Rcurrent, B = 1) %>% drop()
      rs<- Rs[lower.tri(Rs)]
      
      Psi_Rs<- metaSEM::asyCov(Rs, n= asy.n)
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
  reformat<- function(matrixlist){
    k <- ncol(matrixlist[[1]])
    mat_strings <- sapply(matrixlist, function(x) paste(x, collapse = ","))
    freq <- table(mat_strings)
    prop.table(freq)
    
    # Recover matrix
    unique_mats <- lapply(names(freq), function(s) {
      matrix(as.numeric(strsplit(s, ",")[[1]]), nrow = k)
    })
    
    result <- data.frame(
      matrix = I(unique_mats),
      count = as.vector(freq),
      prob = as.vector(prop.table(freq))
    ) %>% as.data.table
    
    result<- result[order(-prob),]
    return(result)
  }
  
}



{
  # Prep: Generate trueR
  p=3L; ad=2; N= 50000L; NCHAIN = 1e4
  hyper_alpha= 0.01
  
  df<- m<- p*(p-1)/2  # degree of freedom
  Target <- er_dag_py(p = p, ad = ad,)
  G0 <- Target$G
  R0 <- Target$R ; kappa(R0)
  B<- Target$B; O<- Target$O
}
{
  asy.n= 200L
  X <- data_gen_py(B, O, n=asy.n, K= 1L)[1,,]
  S <- cor(X); kappa(S)
  if((p*(p-1)/2) > asy.n) warning("high-dimension: p > n")
}
{
  # Truth
  G0
  # Unobtainable Ideal
  Ghat.pc<- pcalg::pc(suffStat= list(C = S, n = asy.n), 
                      indepTest = gaussCItest, 
                      p= p, alpha = as.numeric(hyper_alpha)) %>%
    as(.,"matrix") %>%
    t()
  Ghat.boss<- boss_py(S, asy.n, discount = 1)

}
require(Rgraphviz)
plot(as(G0, "graphNEL"), main = "True DAG")
dev.off()


{
  # mcmc samples
  chain.p3<- make_mcmc(S=S, asy.n= asy.n, NCHAIN=NCHAIN, acceptance.target = 0.4, M= 100)
  
  # bootstrap samples
  bootR <- lapply(1:NCHAIN, function(i){
    boot.id <- sample(asy.n, asy.n, replace = T)
    cor(X[boot.id,])
  })
  bootR<- simplify2array(bootR)
  
  # Learn CDA
  
  fit1_pc<- lapply(1:NCHAIN, function(i){
    fit_obj<- pcalg::pc(suffStat= list(C = bootR[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = as.numeric(hyper_alpha)) 
    est_adj<- t(as(fit_obj, "matrix"))
    return(est_adj)
  })
  fit2_pc<- lapply(1:NCHAIN, function(i){
    fit_obj<- pcalg::pc(suffStat= list(C = chain.p3[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = as.numeric(hyper_alpha)) 
    est_adj<- t(as(fit_obj, "matrix"))
    return(est_adj)
  })
  
  fitboot_boss<- lapply(1:NCHAIN, function(i){boss_py(bootR[,,i], asy.n, discount = 1)})
  fitboot_boss_cpdag <- lapply(1:NCHAIN, function(i){boss_cpdag_py(bootR[,,i], asy.n, discount = 1)})
  fit2_boss<- lapply(1:NCHAIN, function(i){boss_py(chain.p3[,,i], asy.n, discount = 1)})                     
}
dev.off()

{
  
  out1_pc<- reformat(fit1_pc)
  out2_pc<- reformat(fit2_pc)
  out1_boss<- reformat(fit1_boss)
  out2_boss<- reformat(fit2_boss)
  
  k<- max(nrow(out1_pc), nrow(out2_pc), nrow(out1_boss), nrow(out2_boss))
  #k=5
  q1<- k-nrow(out1_pc); q2<- k-nrow(out2_pc); q3<- k- nrow(out1_boss); q4<- k- nrow(out2_boss)
  
  dev.off()
  par(oma = c(0, 0, 6, 0)) 
  par(mfrow=c(1,k))
  
  require(Rgraphviz)
  plot(as(G0, "graphNEL"), main = "True DAG")
  if(k>1) for(i in 2:k) plot.new()
  
  for(i in 1:min(k,nrow(out1_pc)))  plot(as(out1_pc$matrix[[i]], "graphNEL"), main = paste0("freq=",out1_pc[i,2]))
  if(q1 > 0) for(i in (nrow(out1_pc)+1):k) plot.new()
  
  for(i in 1:min(k,nrow(out2_pc)))  plot(as(out2_pc$matrix[[i]], "graphNEL"), main = paste0("freq=",out2_pc[i,2]))
  if(q2 > 0) for(i in (nrow(out2_pc)+1):k) plot.new()
  
  for(i in 1:min(k,nrow(out1_boss)))  plot(as(out1_boss$matrix[[i]], "graphNEL"), main = paste0("freq=",out1_boss[i,2]))
  if(q3 > 0) for(i in (nrow(out1_boss)+1):k) plot.new()
  
  for(i in 1:min(k,nrow(out2_boss)))  plot(as(out2_boss$matrix[[i]], "graphNEL"), main = paste0("freq=",out2_boss[i,2]))
  if(q4 > 0) for(i in (nrow(out2_boss)+1):k) plot.new()
  
  mtext(paste0("Sample size=",asy.n), side = 3, outer = TRUE, line = 4, cex = 1.4, font = 2)
  
}


dev.off()

ncombinations<- sapply(list(out1_pc, out2_pc, out1_boss, out2_boss), nrow)
top5_cont<- sapply(list(out1_pc, out2_pc, out1_boss, out2_boss), function(G) sum(G[1:5,prob]))
summary<- data.frame(ncombinations, top5_cont)
row.names(summary)<- c("boot_pc", "mcmc_pc", "boot_boss", "mcmc_boss")
summary




########################################
# SCALE UP
p <- c(5L)
ad <- 2L
sample.size <- c(50, 100, 200, 500, 1000, 2000, 5000)
ITER <- 10000
hyper_alpha = 0.01
hyper_discount = 1

start = Sys.time()
start
# 1: Generate parameters
{ 
  param0 <- data.table(expand.grid(p = p, ad= ad))
  param0[, out:= mapply(
    function(p, ad) er_dag_py(p = p, ad = ad),
    p, ad,
    SIMPLIFY = FALSE
  )]
  param0[, `:=`(
    G = lapply(out, `[[`, "G"),
    R = lapply(out, `[[`, "R"),
    B = lapply(out, `[[`, "B"),
    O = lapply(out, `[[`, "O")
  )]
  param0[, out:= NULL]
  param0[, graphid:= 1:.N]
}

#2: Resample
{ # parallel run
  #plan(sequential)
  #plan(multisession, workers = min(dim(param)[1], parallelly::availableCores() - 2))
  
  ### the "unattainable ideal": X_
  ### for each graph, generate ITER set of data using the largest sample size
  ### Within the same graph but smaller sample size, subsample (without replacement) from the largest set of sample in each iteration
  param0[, X_:= mapply( 
    function(B,O) data_gen_py(B, O, n= as.integer(max(sample.size)), K=as.integer(ITER)),
    B, O, 
    SIMPLIFY = FALSE
  )]

  param <- param0[rep(1:.N, each = length(sample.size))]
  param[, rowid := .I]
  param[, ss := rep(sample.size, .N / length(sample.size))]
  param[, df:= p*(p-1)/2]
  param<- param[!(ss < df),]  # exclude high dimensional case
  
 
  subsample_array <- function(X, ss) {
    ITER <- dim(X)[1]
    n    <- dim(X)[2]
    p    <- dim(X)[3]
    
    out <- array(NA_real_, dim = c(ITER, ss, p))
    
    for (i in seq_len(ITER)) {
      id <- sample.int(n, ss, replace = FALSE)
      out[i, , ] <- X[i, id, , drop = FALSE]
    }
    
    out
  }
  param[, X_sub := lapply(
    seq_len(.N),
    function(j) subsample_array(X_[[j]], ss[[j]])
  )]
  
  param<- param[,.(graphid, rowid, ss, X_sub)]

  # Take first iteration of X_sub as the "real observed data"
  # Bootstrap will be base don this "real observed data"
  # Store correlation matrices for efficiency
  # We will have three types of sample: A) unobtainable ideal (X_sub), B) bootstrap,  C) mcmc
  param[, S:= lapply(X_sub, function(x) cor(x[1,,]))]
}
{
  ### unobtainable ideal
  ideal_dt <- param[, .(graphid, rowid, X= X_sub)]
  ideal_dt[, `:=`(
    type = "ideal",
    Rs = lapply(X, function(X) {
      n_iter <- dim(X)[1]
      p      <- dim(X)[3]
      out    <- array(NA_real_, c(p, p, n_iter))
      
      for (i in seq_len(n_iter)) {
        out[, ,i] <- cor(X[i, , ])
      }
      out
    })
    )]  
  
  ### bootstrap samples centered around S
  
  boot_dt <- param[, .(graphid, rowid, X = lapply(X_sub, function(x) x[1,,]))]
  boot_dt[, `:=`(
    type = "boot",
    Rs = future_lapply(X, make_boot, ITER = ITER, future.seed= TRUE)
  )]

  ### MCMC samples, starting from a single S
  mcmc_dt <- param[, .(graphid, rowid, S, ss)]
  mcmc_dt[, `:=`(
    type = "mcmc",
    Rs = future_mapply(
      function(S, asy.n) make_mcmc(S, asy.n, ITER),
      S, ss,
      SIMPLIFY = FALSE,
      future.seed = TRUE
    )
  )]
  
  resamples <- rbindlist(list(
    ideal_dt[,.(graphid, rowid, type, Rs)],
    boot_dt[, .(graphid, rowid, type, Rs)],
    mcmc_dt[, .(graphid, rowid, type, Rs)]
  ), use.names = TRUE)
  
  resamples[, ss:= rep(sample.size, length= nrow(resamples))]
  
  # param <- param[resamples, on = "rowid"]
}

print(object.size(resamples), units="auto")
print(object.size(param), units="auto")


{ # sequential run
  param[, rowid := .I]
  resamples <- rbindlist(list(
    param[, .(
      rowid,
      type = "boot",
      Rs = lapply(X, make_boot, ITER = ITER)
    )],
    param[, .(
      rowid,
      type = "mcmc",
      Rs = mapply(
        function(S, asy.n) (make_mcmc(S, asy.n, ITER)),
        S, ss,
        SIMPLIFY = FALSE
      )
    )]
  ), use.names = TRUE)
  param <- param[resamples, on = "rowid"]
}

plan(sequential)
dim(param)

#3: CDA
{ # parallel run
  #options(future.globals.onReference = "ignore")
  plan(multisession, workers = min(dim(resamples)[1], parallelly::availableCores() - 2))
  
  res_list <- future_lapply(seq_len(nrow(resamples)), function(i) {
    library(data.table)
    library(pcalg)
    library(reticulate)
    # initialize Python inside THIS worker
    use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python",required = TRUE)
    reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
    reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
    
    row  <- resamples[i]
    Rs_i <- row$Rs[[1]]
    ss_i <- row$ss
    B    <- dim(Rs_i)[3]
    p_i  <- dim(Rs_i)[1]
    
    boss_out <- lapply(seq_len(B), function(b) {
      ans <- boss_py(
        Rs_i[, , b],
        as.integer(ss_i),
        discount = as.numeric(hyper_discount)
      )
      ans <- reticulate::py_to_r(ans)
      unserialize(serialize(ans, NULL))
    })
    pc_out <- lapply(seq_len(B), function(b) {
      fit_obj <- pcalg::pc(
        suffStat  = list(C = Rs_i[, ,b], n = ss_i/2), ##### ess 
        indepTest = gaussCItest,
        p         = p_i,
        alpha     = as.numeric(hyper_alpha)
      )
      as(fit_obj, "matrix")
    })
    
    rbindlist(list(
      data.table(
        graphid = row$graphid,
        rowid   = row$rowid,
        type    = row$type,
        ss    = ss_i,
        method  = "boss",
        out     = list(boss_out)
      ),
      data.table(
        graphid = row$graphid,
        rowid   = row$rowid,
        type    = row$type,
        ss    = ss_i,
        method  = "pc",
        out     = list(pc_out)
      )
    ), use.names = TRUE)
  }, future.seed = TRUE)
  
  res <- rbindlist(res_list, use.names = TRUE)
  plan(sequential)
}
print(object.size(res), units= "auto")


{ # sequential run
  boss_dt <- param[, {
    Rs_i <- Rs[[1]]
    ss_i <- ss
    .(
      method = "boss",
      out = list(
        lapply(seq_len(ITER), function(b) {
          boss_py(Rs_i[, , b], as.integer(ss_i), discount = as.numeric(hyper_discount))
          })
      )
    )
  }, by = .(graphid, rowid, p, ad, ss, type)]
  
  pc_dt <- param[, {
    Rs_i <- Rs[[1]]
    ss_i <- ss
    p <- dim(Rs_i)[1]
    .(
      method = "pc",
      out = list(
        lapply(seq_len(ITER), function(b) {
          fit_obj <- pcalg::pc(suffStat = list(C = Rs_i[,,b], n = ss_i), indepTest = gaussCItest, p = p, alpha = as.numeric(hyper_alpha))
          (as(fit_obj, "matrix"))        
          })
      )
    )
  }, by = .(graphid, rowid, p, ad, ss, type)]
  
  res <- rbindlist(list(boss_dt, pc_dt), use.names = TRUE)
}
gc()


# convert into cpdag
{
  res<- res[param0[, .(G = list(G[[1]]), p, ad), by = graphid], on=.(graphid)]
  res_cpdag <- copy(res)
  res_cpdag[, out_cpdag := lapply(out, function(X) {
    lapply(X, function(dag) 1L*dag2cpdag(dag))
  })]
  res_cpdag[, out:=NULL]
  setnames(res_cpdag, "out_cpdag", "out")
  #res_cpdag
}


#4: Enumerate combination
{
  res[,reformat_out:= lapply(out, reformat) ]
  res[,num.graphs:= lapply(reformat_out, nrow)]
  res[,top5_cont:= sapply(reformat_out, function(X){
    id= min(5, nrow(X))
    sum(X[1:id,prob])
  })]
  
  res_cpdag[,reformat_out:= lapply(out, reformat) ]
  res_cpdag[,num.graphs:= lapply(reformat_out, nrow)]
  res_cpdag[,top5_cont:= sapply(reformat_out, function(X){
    id= min(5, nrow(X))
    sum(X[1:id,prob])
  })]
}

res[, out2:= out]
res[method=="boss", out2:= lapply(out, function(X) {
  lapply(X, function(x) 1L*dag2cpdag(x))
})]


res_p3<- res
res_cpdag_p3 <- res_cpdag

dcast.data.table(data= res, formula = method + graphid + rowid + p + ad + ss  ~  type, value.var = "num.graphs" )
dcast.data.table(data= res_cpdag, formula = method  +graphid + rowid + p + ad + ss ~  type, value.var = "num.graphs" )









#5: Power Analysis
### Add in Ghat: the estimated G from S using the corresponding method
res1<- res[type=="ideal", .(Ghat= out[[1]][1]), by=.(graphid, rowid, method)]
res<- res[res1, on=.(graphid, rowid, method)]

res_cpdag<- res_cpdag[res1[, Ghat:=lapply(Ghat, function(dag) 1L*dag2cpdag(dag))], on=.(graphid, rowid, method)]


{ # parallel
  plan(multisession, workers = parallelly::availableCores() - 1)
  metric_list <- future_Map(function(out_i, G_i) {
    perf_i <- sapply(out_i, function(K) {
      tmp <- eval_cda(G_i, K)
      c(tmp$conf.adj, tmp$conf.orient)
    }) |> t()
    
    list(
      adj_precision   = perf_i[, 1] / (perf_i[, 1] + perf_i[, 3]),
      adj_sensitivity = perf_i[, 1] / (perf_i[, 1] + perf_i[, 2]),
      ort_precision   = perf_i[, 5] / (perf_i[, 5] + perf_i[, 7]),
      ort_sensitivity = perf_i[, 5] / (perf_i[, 5] + perf_i[, 6])
    )
  }, res$out, res$G, future.seed = TRUE)
  
  
  metric_list_cpdag <- future_Map(function(out_i, G_i) {
    perf_i <- sapply(out_i, function(K) {
      tmp <- eval_cda(G_i, K)
      c(tmp$conf.adj, tmp$conf.orient)
    }) |> t()
    
    list(
      adj_precision   = perf_i[, 1] / (perf_i[, 1] + perf_i[, 3]),
      adj_sensitivity = perf_i[, 1] / (perf_i[, 1] + perf_i[, 2]),
      ort_precision   = perf_i[, 5] / (perf_i[, 5] + perf_i[, 7]),
      ort_sensitivity = perf_i[, 5] / (perf_i[, 5] + perf_i[, 6])
    )
  }, res_cpdag$out, res_cpdag$G, future.seed = TRUE)
  
  res[, `:=`(
    adj_precision   = lapply(metric_list, `[[`, "adj_precision"),
    adj_sensitivity = lapply(metric_list, `[[`, "adj_sensitivity"),
    ort_precision   = lapply(metric_list, `[[`, "ort_precision"),
    ort_sensitivity = lapply(metric_list, `[[`, "ort_sensitivity")
  )]
  
  res_cpdag[, `:=`(
    adj_precision   = lapply(metric_list_cpdag, `[[`, "adj_precision"),
    adj_sensitivity = lapply(metric_list_cpdag, `[[`, "adj_sensitivity"),
    ort_precision   = lapply(metric_list_cpdag, `[[`, "ort_precision"),
    ort_sensitivity = lapply(metric_list_cpdag, `[[`, "ort_sensitivity")
  )]
  plan(sequential)
}


{ #sequential
  res[, perf := future_Map(function(out_i, G_i) {
    perf_i <- sapply(out_i, function(K) {
      tmp <- eval_cda(G_i, K)
      c(tmp$conf.adj, tmp$conf.orient)
    }) |> t()
    
    colnames(perf_i) <- c(
      paste0("adj_", c("TP", "FN", "FP", "TN")),
      paste0("ort_", c("TP", "FN", "FP", "TN"))
    )
    perf_i
  }, out, G)]
  
  res[,`:=`(
    adj_precision= lapply(perf, function(X)  X[,1]/(X[,1]+X[,3])),
    adj_sensitivity= lapply(perf, function(X)  X[,1]/(X[,1]+X[,2])), 
    ort_precision = lapply(perf, function(X)  X[,5]/(X[,5]+X[,7])),
    ort_sensitivity = lapply(perf, function(X)  X[,5]/(X[,5]+X[,6]))
  )]
}
end= Sys.time()


gc(); getwd()


saveRDS(list(res = res, param = param), "UseCase/run1_p9to10.rds")
tmp <- readRDS("UseCase/run1_p3to8.rds")
res   <- tmp$res
param <- tmp$param

##############

plan(sequential)
res[, .(
  adj_precision_med   = mean(unlist(adj_precision), na.rm = TRUE),
  adj_sensitivity_med = mean(unlist(adj_sensitivity), na.rm = TRUE),
  ort_precision_med   = mean(unlist(ort_precision), na.rm = TRUE),
  ort_sensitivity_med = mean(unlist(ort_sensitivity), na.rm = TRUE)
), by = .(type, method, p, ad, ss)][ss==1000,]
res_cpdag[, .(
  adj_precision_med   = mean(unlist(adj_precision), na.rm = TRUE),
  adj_sensitivity_med = mean(unlist(adj_sensitivity), na.rm = TRUE),
  ort_precision_med   = mean(unlist(ort_precision), na.rm = TRUE),
  ort_sensitivity_med = mean(unlist(ort_sensitivity), na.rm = TRUE)
), by = .(type, method, p, ad, ss)][ss==1000 & method=="boss",]

