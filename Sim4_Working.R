# PerturbR Sim3
rm(list=ls()); gc()
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr)
use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
py_config() 

source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R")  # Wald cloud
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R")  # Base CDA Functions

# Source python files
source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")



# Parameters --------------------------------------------------------------

### Sim parameters
# p : 10, 20, 30. asy.n: 500
{
  nTarget <- 200L
  p <- 20L
  ad <- 2L
  asy.n <- 300L
  ITER <- 300L
  alpha2 <- 0.01
  main_title <- paste0(
    nTarget, " Graphs (p=", p, ", ad=", ad, 
    "), each with ", ITER, " resamples, n=", asy.n, 
    ", alpha=", alpha2
  )
  out_dir<- "~/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/plots/Performance/"
}

### Default Functions
{
  methods <- c("waldcloud", "unifcloud_001", "unifcloud_01", "unifcloud_02", "bootstrap")
  n_obs_by_method <- c(
    waldcloud = asy.n,
    unifcloud_001 = asy.n,
    unifcloud_01 = asy.n,
    unifcloud_02 = asy.n,
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
  
  
  run_pc_mat <- function(Cmat, n_obs, alpha) {
    p= ncol(Cmat)
    as(
      pcalg::pc(
        suffStat = list(C = Cmat, n = n_obs),
        indepTest = gaussCItest,
        p = p,
        alpha = alpha
      ),
      "matrix"
    )
  }
  
  get_metric <- function(G0, Ghat) {
    out <- eval_cda(true_adj = t(G0), est_adj = t(Ghat))
    list(adj = out$metric.adj, orient = out$metric.orient)
  }
  
}




# Simulation --------------------------------------------------------------

# Parallel Run ------------------------------------------------------------
{
  {
    library(parallelly); library(future.apply)
    plan(multisession, workers = max(1, parallelly::availableCores() - 1))
    run_one_nT <- function(nT, p, ad, asy.n, ITER, methods, n_obs_by_method, alpha2) {
      # Worker-local init
      reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
      reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
      source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
      source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)
      
      Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 1L)
      G0 <- Target$G
      R0 <- Target$R
      X <- Target$X[1, , ]
      Rhat <- cor(X)
      
      # keep same top-level structure
      out <- list(
        truth = G0,
        optimal = list(
          pc = run_pc_mat(R0, asy.n, alpha = alpha2),
          boss = boss_py(R0, asy.n)
        ),
        baseline = list(
          pc = run_pc_mat(Rhat, asy.n, alpha = alpha2),
          boss = boss_py(Rhat, asy.n)
        ),
        methods = setNames(vector("list", length(methods)), methods)
      )
      
      # build Rhat_arr method-by-method with error capture
      rhat_builders <- list(
        waldcloud     = function() rwaldcloud(Rhat, n = asy.n, B = ITER),
        unifcloud_001 = function() runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.01),
        unifcloud_01  = function() runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.1),
        unifcloud_02  = function() runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.2),
        unifcloud_095  = function() runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.95),
        bootstrap     = function() boot_cor(X, asy.n = asy.n, B = ITER)
      )
      
      Rhat_arr <- setNames(vector("list", length(methods)), methods)
      Rhat_err <- setNames(vector("list", length(methods)), methods)
      
      for (m in methods) {
        tmp <- tryCatch(rhat_builders[[m]](), error = function(e) e)
        if (inherits(tmp, "error")) {
          Rhat_arr[[m]] <- NULL
          Rhat_err[[m]] <- conditionMessage(tmp)
        } else {
          Rhat_arr[[m]] <- tmp
          Rhat_err[[m]] <- NULL
        }
      }
      
      # fill out$methods with same nesting; skip fit calls for failed methods
      for (m in methods) {
        if (!is.null(Rhat_err[[m]])) {
          out$methods[[m]] <- list(
            pc = NULL,
            boss = NULL,
            error = Rhat_err[[m]]
          )
          next
        }
        
        n_obs <- n_obs_by_method[[m]]
        out$methods[[m]] <- list(
          pc = lapply(seq_len(ITER), function(i) run_pc_mat(Rhat_arr[[m]][, , i], n_obs, alpha = alpha2)),
          boss = lapply(seq_len(ITER), function(i) boss_py(Rhat_arr[[m]][, , i], n_obs)),
          error = NULL
        )
      }
      
      out
    }
    
    res <- future_lapply(
      seq_len(nTarget),
      run_one_nT,
      p = p, ad = ad, asy.n = asy.n, ITER = ITER, alpha2= alpha2,
      methods = methods, n_obs_by_method = n_obs_by_method,
      future.packages = c("pcalg", "reticulate"),
      future.seed = TRUE,
      future.globals = c("run_pc_mat")  # keep globals minimal
    )
    
    plan(sequential)
    
  }
  {
    rebuild_fit_from_res <- function(res, drop_error_nT = TRUE) {
      if (length(res) == 0L) stop("`res` is empty.")
      
      fixed_keys <- c("truth", "optimal", "baseline", "methods")
      
      # Detect if an nT has any method$error != NULL
      has_method_error <- function(x) {
        method_block <- if ("methods" %in% names(x) && is.list(x$methods)) {
          x$methods
        } else {
          x[setdiff(names(x), fixed_keys)]
        }
        
        if (length(method_block) == 0L) return(FALSE)
        
        any(vapply(
          method_block,
          function(mm) is.list(mm) && !is.null(mm$error),
          logical(1)
        ))
      }
      
      dropped_nT <- integer(0)
      kept_nT <- seq_along(res)
      
      if (drop_error_nT) {
        dropped_nT <- which(vapply(res, has_method_error, logical(1)))
        kept_nT <- setdiff(seq_along(res), dropped_nT)
        res <- res[kept_nT]
      }
      
      nTarget <- length(res)
      if (nTarget == 0L) stop("All nT were dropped due to errors.")
      
      # infer methods from first kept element
      methods <- if ("methods" %in% names(res[[1]]) && is.list(res[[1]]$methods)) {
        names(res[[1]]$methods)
      } else {
        setdiff(names(res[[1]]), fixed_keys)
      }
      
      fit <- list(
        truth = vector("list", nTarget),
        optimal = list(
          pc = vector("list", nTarget),
          boss = vector("list", nTarget)
        ),
        baseline = list(
          pc = vector("list", nTarget),
          boss = vector("list", nTarget)
        )
      )
      
      for (m in methods) fit[[m]] <- vector("list", nTarget)
      
      for (i in seq_len(nTarget)) {
        x <- res[[i]]
        
        fit$truth[[i]] <- x$truth
        fit$optimal$pc[[i]] <- x$optimal$pc
        fit$optimal$boss[[i]] <- x$optimal$boss
        fit$baseline$pc[[i]] <- x$baseline$pc
        fit$baseline$boss[[i]] <- x$baseline$boss
        
        for (m in methods) {
          fit[[m]][[i]] <- if ("methods" %in% names(x) && is.list(x$methods)) x$methods[[m]] else x[[m]]
        }
      }
      
      attr(fit, "dropped_nT") <- dropped_nT
      attr(fit, "kept_nT_original_index") <- kept_nT
      fit
    }
    fit1 <- rebuild_fit_from_res(res)
  }
  { # Error Count
    # nT indices that have at least one method error
    nT_with_error <- which(vapply(
      res,
      function(x) any(vapply(x$methods, function(mm) !is.null(mm$error), logical(1))),
      logical(1)
    ))
    nT_with_error
  }
  
  # Update nT: 
  nTarget = nTarget- length(nT_with_error)
}



#fit$optimal: learn based on R0: G0hat
#fit$baseline: learn based on Rhat
#fit$waldcloud and other methods.... learn based on perturbation of Rhat

# ETL and Plotting --------------------------------------------------------
### Table ETL 
{# New
  library(data.table)
  
  safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
  build_eval_table <- function(
    fit, ITER, methods, algos, targets, nTarget,
    parallel = FALSE, n_workers = max(1L, parallel::detectCores(logical = FALSE) - 1L)
  ) {
    # local alias to avoid repeated global lookup
    e_cda <- eval_cda
    
    # one work unit = (algo, nT)
    jobs <- data.table::CJ(
      algo = algos,
      nT = seq_len(nTarget),
      sorted = FALSE
    )

    eval_one_job <- function(algo, nT) {
      G0 <- fit$truth[[nT]]
      G0hat <- fit$optimal[[algo]][[nT]]
      Ghat_base <- fit$baseline[[algo]][[nT]]
      true_map <- list(G0 = G0, G0hat = G0hat, Ghat = Ghat_base)
      
      # pre-size roughly to reduce reallocations
      approx_n <- length(targets) * (sum(methods == "baseline") + ITER * sum(methods != "baseline"))
      out <- vector("list", approx_n)
      r <- 0L
      
      for (m in methods) {
        if (m == "baseline") {
          # baseline uses one estimate only
          G_est <- Ghat_base
          if (!is.null(G_est)) {
            for (tgt in targets) {
              G_true <- true_map[[tgt]]
              out_mat <- e_cda(true_adj = G_true, est_adj = G_est)$out
              r <- r + 1L
              out[[r]] <- list(
                method = m, algo = algo, target = tgt, nT = nT, i = 1L,
                TP_adj = out_mat[1, 2], FP_adj = out_mat[2, 2], FN_adj = out_mat[3, 2], TN_adj = out_mat[4, 2],
                TP_ort = out_mat[1, 3], FP_ort = out_mat[2, 3], FN_ort = out_mat[3, 3], TN_ort = out_mat[4, 3]
              )
            }
          }
        } else {
          fit_m_nT_algo <- fit[[m]][[nT]][[algo]]
          if (!is.null(fit_m_nT_algo)) {
            n_i <- min(ITER, length(fit_m_nT_algo))
            for (i in seq_len(n_i)) {
              G_est <- fit_m_nT_algo[[i]]
              if (is.null(G_est)) next
              
              for (tgt in targets) {
                G_true <- true_map[[tgt]]
                out_mat <- e_cda(true_adj = G_true, est_adj = G_est)$out
                r <- r + 1L
                out[[r]] <- list(
                  method = m, algo = algo, target = tgt, nT = nT, i = i,
                  TP_adj = out_mat[1, 2], FP_adj = out_mat[2, 2], FN_adj = out_mat[3, 2], TN_adj = out_mat[4, 2],
                  TP_ort = out_mat[1, 3], FP_ort = out_mat[2, 3], FN_ort = out_mat[3, 3], TN_ort = out_mat[4, 3]
                )
              }
            }
          }
        }
      }
      
      if (r == 0L) return(NULL)
      data.table::rbindlist(out[seq_len(r)], use.names = TRUE, fill = FALSE)
    }
    
    #for (j in 1:nrow(jobs)){
    #  eval_one_job(jobs$algo[j], jobs$nT[j])
    #}


    if (!parallel) {
      chunks <- lapply(seq_len(nrow(jobs)), function(j) {
        eval_one_job(jobs$algo[j], jobs$nT[j])
      })
    } else {
      # parallel over independent (algo, nT) jobs
      chunks <- parallel::mclapply(
        X = seq_len(nrow(jobs)),
        FUN = function(j) eval_one_job(jobs$algo[j], jobs$nT[j]),
        mc.cores = n_workers
      )
    }
    
    dt <- data.table::rbindlist(chunks, use.names = TRUE, fill = TRUE)
    
    dt[, `:=`(
      precision_adj = safe_div(TP_adj, TP_adj + FP_adj),
      recall_adj    = safe_div(TP_adj, TP_adj + FN_adj),
      precision_ort = safe_div(TP_ort, TP_ort + FP_ort),
      recall_ort    = safe_div(TP_ort, TP_ort + FN_ort)
    )]
    dt[, `:=`(
      f1_adj = safe_div(2 * precision_adj * recall_adj, precision_adj + recall_adj),
      f1_ort = safe_div(2 * precision_ort * recall_ort, precision_ort + recall_ort)
    )]
    
    return(dt)
  }
  
  dt_eval <- build_eval_table(
    fit = fit1,
    ITER= ITER,
    methods = c("baseline","bootstrap","waldcloud","unifcloud_001","unifcloud_01","unifcloud_02"),
    algos = c("pc","boss"),
    targets = c("G0","G0hat"),
    nTarget = nTarget, 
    parallel = TRUE, 
    n_workers = 8
  )

}


## Plotting
### Functions
{ # Plot
  
  transform_dt<- function(df_eval){
    dt_plot <- melt(
      df_eval,
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
    
    return(dt_nt_median)
    
  }
  
  plot_nt_median_box <- function(dt_nt_median,
                                 algo_pick = "pc",
                                 eval_type_pick = "adj",
                                 metric_pick = "precision",
                                 method_levels,
                                 x_tick,
                                 target_levels= c("G0", "G0hat")) {
    
    df <- dt_nt_median[
      algo == algo_pick &
        eval_type == eval_type_pick &
        metric == metric_pick
    ]
    
    df[, method := factor(method, levels = method_levels)]
    df[, target := factor(target, levels = c("G0", "G0hat"), labels = c("Compared to G0", "Compared to G0hat"))]
    
    ggplot(df, aes(x = method, y = score_med, fill = target)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.7, show.legend = TRUE) +
      #facet_wrap(~ target, nrow = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      #coord_flip() +
      labs(
        x = "Method",
        y = NULL,
        #y = paste0("Median ", metric_pick, " (per graph)"),
        title = paste0(toupper(algo_pick), " — ", eval_type_pick, " ", metric_pick, " (Median per graph)"),
        fill = "Compared to"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      scale_x_discrete(limits = method_levels, labels = parse(text = x_tick))

  }
  save_png_with_index <- function(plot,out_dir,main_title,prefix = "BOSS") {
    # Ensure directory exists
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    # List existing PNG files
    existing <- list.files(out_dir, pattern = "\\.png$", full.names = FALSE)
    
    # Extract leading numbers (only if at start of filename)
    nums <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", existing)))
    nums <- nums[!is.na(nums)]
    
    # Compute next index
    next_index <- if (length(nums) == 0) 0 else max(nums) + 1
    
    # Clean title
    clean_title <- gsub("[^[:alnum:]]+", "", main_title)
    
    # Build filename
    filename <- file.path(
      out_dir,
      paste0(next_index, "_", prefix, " (", clean_title, ").png")
    )
    
    # Save plot
    ggplot2::ggsave(filename = filename, plot = plot)
    
    invisible(filename)
  }
  
}

### Draw and Save
{
  dt_nt_median<- transform_dt(dt_eval)
  dt_nt_median<- dt_nt_median[target!="Ghat"]
  methods.level<- c("baseline", "bootstrap", "waldcloud", "unifcloud_001", "unifcloud_01", "unifcloud_02")
  #methods.display <- c("noperturb", "boot","$W(\\hat{R})$","$U_{0.99}(\\hat{R})$","$U_{0.9}(\\hat{R})$","$U_{0.8}(\\hat{R})$")
  methods.display <- c("no~perturb","boot","W(hat(R))","U[0.99](hat(R))","U[0.9](hat(R))","U[0.8](hat(R))")
  
  p1<- plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="adj", metric_pick="precision", methods.level, x_tick = methods.display) 
  p2<- plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="adj", metric_pick="recall", methods.level, x_tick = methods.display)
  p3<- plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="adj", metric_pick="f1", methods.level, x_tick = methods.display)
  p4<- plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="orient", metric_pick="precision", methods.level, x_tick = methods.display)
  p5<- plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="orient", metric_pick="recall", methods.level, x_tick = methods.display)
  p6<- plot_nt_median_box(dt_nt_median, algo_pick="pc", eval_type_pick="orient", metric_pick="f1",methods.level, x_tick = methods.display)
  fig <- ggarrange(p1, p4,  p2, p5,  p3, p6,  ncol = 2, nrow = 3,  align = "h", common.legend = TRUE,  legend = "top" )
  fig <- ggpubr::annotate_figure(  fig,  top = text_grob(main_title,face = "bold", size = 14))
  fig
  save_png_with_index(plot = fig,out_dir = out_dir,main_title = main_title, prefix="PC Individual"  )
  
  p1<- plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="adj", metric_pick="precision",methods.level, x_tick = methods.display)
  p2<- plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="adj", metric_pick="recall",methods.level, x_tick = methods.display)
  p3<- plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="adj", metric_pick="f1",methods.level, x_tick = methods.display)
  p4<- plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="orient", metric_pick="precision",methods.level, x_tick = methods.display)
  p5<- plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="orient", metric_pick="recall",methods.level, x_tick = methods.display)
  p6<- plot_nt_median_box(dt_nt_median, algo_pick="boss", eval_type_pick="orient", metric_pick="f1",methods.level, x_tick = methods.display)
  #ggpubr::ggarrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2)
  
  fig <- ggarrange(p1, p4,  p2, p5,  p3, p6,  ncol = 2, nrow = 3,  align = "h", common.legend = TRUE,  legend = "top" )
  fig <- ggpubr::annotate_figure(  fig,  top = text_grob(main_title,face = "bold", size = 14))
  fig
  save_png_with_index(plot = fig,out_dir = out_dir,main_title = main_title, prefix="BOSS Individual"  )
  
}



# ENSEMBLE ----------------------------------------------------------------

# Ensemble results
threshold= 0.6  # stability > threshold , then convert edge to 1
{
  methods = c("bootstrap","waldcloud","unifcloud_001","unifcloud_01","unifcloud_02")
  fit.ensemble= NULL
  fit.ensemble$truth<- fit1$truth
  fit.ensemble$optimal<- fit1$optimal
  fit.ensemble$baseline<- fit1$baseline
  fit.ensemble[[m]] <- setNames(vector("list", length(methods)), methods)
  for (m in methods) {
    fit.ensemble[[m]] <- vector("list", nTarget)
    
    for (nT in seq_len(nTarget)) {
      pc_list <- fit1[[m]][[nT]][["pc"]]
      boss_list <- fit1[[m]][[nT]][["boss"]]
      
      fit.ensemble[[m]][[nT]] <- list(
        pc   = list(ifelse(Reduce(`+`, pc_list)/length(pc_list)>threshold, 1, 0)),
        boss = list(ifelse(Reduce(`+`, boss_list)/length(pc_list)>threshold, 1, 0))
      )
    }
  }
  
}


dt_eval_ens <- build_eval_table(
  fit = fit.ensemble,
  ITER= ITER,
  methods = c("baseline","waldcloud","unifcloud_001","unifcloud_01","unifcloud_02","bootstrap"),
  algos = c("pc","boss"),
  targets = c("G0","G0hat"),
  nTarget = nTarget, 
  parallel = FALSE, 
  n_workers = 8
)


# 
{ 
  # Plotting
  dt_nt_median_ens<- transform_dt(dt_eval_ens)
  dt_nt_median_ens<- dt_nt_median_ens[target!="Ghat"]
  methods.level<- c("baseline", "bootstrap", "waldcloud", "unifcloud_001", "unifcloud_01", "unifcloud_02")
  
  p1<- plot_nt_median_box(dt_nt_median_ens, algo_pick="pc", eval_type_pick="adj", metric_pick="precision", methods.level, x_tick = methods.display)
  p2<- plot_nt_median_box(dt_nt_median_ens, algo_pick="pc", eval_type_pick="adj", metric_pick="recall", methods.level, x_tick = methods.display)
  p3<- plot_nt_median_box(dt_nt_median_ens, algo_pick="pc", eval_type_pick="adj", metric_pick="f1", methods.level, x_tick = methods.display)
  p4<- plot_nt_median_box(dt_nt_median_ens, algo_pick="pc", eval_type_pick="orient", metric_pick="precision", methods.level, x_tick = methods.display)
  p5<- plot_nt_median_box(dt_nt_median_ens, algo_pick="pc", eval_type_pick="orient", metric_pick="recall", methods.level, x_tick = methods.display)
  p6<- plot_nt_median_box(dt_nt_median_ens, algo_pick="pc", eval_type_pick="orient", metric_pick="f1",methods.level, x_tick = methods.display)
  
  fig <- ggarrange(p1, p4,  p2, p5,  p3, p6,  ncol = 2, nrow = 3,  align = "h", common.legend = TRUE,  legend = "top" )
  fig <- ggpubr::annotate_figure(  fig,  top = text_grob(paste0("Ensemble results: ", main_title),face = "bold", size = 14))
  fig
  save_png_with_index(plot = fig,out_dir = out_dir,main_title = main_title, prefix="PC Ensemble")
  
  p1<- plot_nt_median_box(dt_nt_median_ens, algo_pick="boss", eval_type_pick="adj", metric_pick="precision", methods.level, x_tick = methods.display)
  p2<- plot_nt_median_box(dt_nt_median_ens, algo_pick="boss", eval_type_pick="adj", metric_pick="recall", methods.level, x_tick = methods.display)
  p3<- plot_nt_median_box(dt_nt_median_ens, algo_pick="boss", eval_type_pick="adj", metric_pick="f1", methods.level, x_tick = methods.display)
  p4<- plot_nt_median_box(dt_nt_median_ens, algo_pick="boss", eval_type_pick="orient", metric_pick="precision", methods.level, x_tick = methods.display)
  p5<- plot_nt_median_box(dt_nt_median_ens, algo_pick="boss", eval_type_pick="orient", metric_pick="recall", methods.level, x_tick = methods.display)
  p6<- plot_nt_median_box(dt_nt_median_ens, algo_pick="boss", eval_type_pick="orient", metric_pick="f1",methods.level, x_tick = methods.display)
  
  fig <- ggarrange(p1, p4,  p2, p5,  p3, p6,  ncol = 2, nrow = 3,  align = "h", common.legend = TRUE,  legend = "top" )
  fig <- ggpubr::annotate_figure(  fig,  top = text_grob(paste0("Ensemble results: ", main_title),face = "bold", size = 14))
  fig
  save_png_with_index(plot = fig,out_dir = out_dir,main_title = main_title, prefix="BOSS Ensemble")
  
}





#####
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





## OLD

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










# Full Simulation 2 ---------------------------------------------------


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
    unifcloud_02  = runifcloud(Rhat, n = asy.n, B = ITER, alpha = 0.2),
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


