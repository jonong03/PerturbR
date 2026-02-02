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


# On time Simulation -----------------------------------------------------
p= 20L; ad= 3L; asy.n= 20000L; ITER=300L
ps= p*(p-1)/2
#alpha1<- 0.05   # alpha1= 1 minus coveragelevel
alpha2<- 0.05   # used in CI test
{
  # Truth Layer:
  Target = er_dag_py(p=p, ad=ad, n= asy.n, K=1L)  # Consists of K data sets
  G0 <- Target$G
  #print(G0)
  R0 <- Target$R
  X<- Target$X[1,,]   # The only observed data
  
  # Generate Rsamples. Approaches: 1) Bootstrap, 2) 1-alpha uniform wald cloud
  # Returns 3d-Array
  #Rhat.waldcloud<- rwaldcloud(Target$R, n = asy.n, B = ITER)
  #Rhat.corX <- lapply(1:ITER, function(i) cor(Target$X[i,,])) %>% simplify2array
  Rhat.bootstrap <- lapply(1:ITER, function(i){
    BOOT.ID<- sample(asy.n, asy.n, replace= TRUE)
    X.boot<- X[BOOT.ID,]
    return(cor(X.boot))
  }) %>% simplify2array
  Rhat.unifcloud_001<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.01)
  Rhat.unifcloud_005<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.05)
  Rhat.unifcloud_01<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.1)
  Rhat.unifcloud_02<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.2)
  
  # Run CDA
  # PC
  {
    fit0<- pcalg::pc(suffStat= list(C = R0, n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2 )
    Ghat.0 <- as(fit0, "matrix")
    #fit.waldcloud<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.waldcloud[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    #fit.corX<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.corX[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    fit.unifcloud_001<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_001[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    fit.unifcloud_005<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_005[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    fit.unifcloud_01<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_01[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    fit.unifcloud_02<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_02[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    fit.bootstrap<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.bootstrap[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
    
    #Ghat.waldcloud<- lapply(fit.waldcloud, function(Y) as(Y, "matrix"))
    #Ghat.corX<- lapply(fit.corX, function(Y) as(Y, "matrix"))
    Ghat.unifcloud_001<- lapply(fit.unifcloud_001, function(Y) as(Y, "matrix"))
    Ghat.unifcloud_005<- lapply(fit.unifcloud_005, function(Y) as(Y, "matrix"))
    Ghat.unifcloud_01<- lapply(fit.unifcloud_01, function(Y) as(Y, "matrix"))
    Ghat.unifcloud_02<- lapply(fit.unifcloud_02, function(Y) as(Y, "matrix"))
    Ghat.bootstrap<- lapply(fit.bootstrap, function(Y) as(Y, "matrix"))
    
    # Compare Ghats with G0
    #perf.waldcloud<- sapply(Ghat.waldcloud, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric)) 
    #perf.corX<- sapply(Ghat.corX, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric)) 
    perf_adj.unifcloud_001<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
    perf_adj.unifcloud_005<- sapply(Ghat.unifcloud_005, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
    perf_adj.unifcloud_01<- sapply(Ghat.unifcloud_01, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
    perf_adj.unifcloud_02<- sapply(Ghat.unifcloud_02, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
    perf_adj.bootstrap<- sapply(Ghat.bootstrap, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
    
    perf_ort.unifcloud_001<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    perf_ort.unifcloud_005<- sapply(Ghat.unifcloud_005, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    perf_ort.unifcloud_01<- sapply(Ghat.unifcloud_01, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    perf_ort.unifcloud_02<- sapply(Ghat.unifcloud_02, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    perf_ort.bootstrap<- sapply(Ghat.bootstrap, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    
    # For PC plotting
    my_scale<- scale_x_discrete(
      labels = c("bootstrap" = "Bootstrap", "unifcloud_0.01" = "99% Cloud", "unifcloud_0.05" = "95% Cloud", 
                 "unifcloud_0.1" = "90% Cloud", "unifcloud_0.2" = "80% Cloud")
    )
    {
      df_adjacency_accuracy= bind_rows(
        #waldcloud= data.frame(value= perf.waldcloud[1,], group="waldcloud"),
        #corX= data.frame(value= perf.corX[1,], group="corData"),
        unifcloud_001= data.frame(value= perf_adj.unifcloud_001[1,], group="unifcloud_0.01"),
        unifcloud_005= data.frame(value= perf_adj.unifcloud_005[1,], group="unifcloud_0.05"),
        unifcloud_01= data.frame(value= perf_adj.unifcloud_01[1,], group="unifcloud_0.1"),
        unifcloud_02= data.frame(value= perf_adj.unifcloud_02[1,], group="unifcloud_0.2"),
        boostrap= data.frame(value= perf_adj.bootstrap[1,], group="bootstrap")
      )
      df_adjacency_sensitivity= bind_rows(
        #waldcloud= data.frame(value= perf.waldcloud[1,], group="waldcloud"),
        #corX= data.frame(value= perf.corX[1,], group="corData"),
        unifcloud_001= data.frame(value= perf_adj.unifcloud_001[2,], group="unifcloud_0.01"),
        unifcloud_005= data.frame(value= perf_adj.unifcloud_005[2,], group="unifcloud_0.05"),
        unifcloud_01= data.frame(value= perf_adj.unifcloud_01[2,], group="unifcloud_0.1"),
        unifcloud_02= data.frame(value= perf_adj.unifcloud_02[2,], group="unifcloud_0.2"),
        boostrap= data.frame(value= perf_adj.bootstrap[2,], group="bootstrap")
      )
      df_orientation_accuracy= bind_rows(
        #waldcloud= data.frame(value= perf.waldcloud[1,], group="waldcloud"),
        #corX= data.frame(value= perf.corX[1,], group="corData"),
        unifcloud_001= data.frame(value= perf_ort.unifcloud_001[1,], group="unifcloud_0.01"),
        unifcloud_005= data.frame(value= perf_ort.unifcloud_005[1,], group="unifcloud_0.05"),
        unifcloud_01= data.frame(value= perf_ort.unifcloud_01[1,], group="unifcloud_0.1"),
        unifcloud_02= data.frame(value= perf_ort.unifcloud_02[1,], group="unifcloud_0.2"),
        boostrap= data.frame(value= perf_ort.bootstrap[1,], group="bootstrap")
      )
      df_orientation_sensitivity= bind_rows(
        #waldcloud= data.frame(value= perf.waldcloud[1,], group="waldcloud"),
        #corX= data.frame(value= perf.corX[1,], group="corData"),
        unifcloud_001= data.frame(value= perf_ort.unifcloud_001[2,], group="unifcloud_0.01"),
        unifcloud_005= data.frame(value= perf_ort.unifcloud_005[2,], group="unifcloud_0.05"),
        unifcloud_01= data.frame(value= perf_ort.unifcloud_01[2,], group="unifcloud_0.1"),
        unifcloud_02= data.frame(value= perf_ort.unifcloud_02[2,], group="unifcloud_0.2"),
        boostrap= data.frame(value= perf_ort.bootstrap[2,], group="bootstrap")
      )
      
      g_adj_acc<- ggplot2::ggplot(df_adjacency_accuracy, aes(x = group, y = value)) +
        geom_boxplot() +
        ggtitle("PC: Adjacency Accuracy")+
        theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      g_adj_sen<- ggplot2::ggplot(df_adjacency_sensitivity, aes(x = group, y = value)) +
        geom_boxplot() +
        ggtitle("PC: Adjacency Sensitivity")+
        theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      g_ort_acc<- ggplot2::ggplot(df_orientation_accuracy, aes(x = group, y = value)) +
        geom_boxplot() +
        ggtitle("PC: Orientation Accuracy")+
        theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      g_ort_sen<- ggplot2::ggplot(df_orientation_sensitivity, aes(x = group, y = value)) +
        geom_boxplot() +
        ggtitle("PC: Orientation Sensitivity")+
        theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      ggpc<- ggpubr::ggarrange(g_adj_acc, g_adj_sen,g_ort_acc, g_ort_sen,ncol = 2, nrow=2)
      
    }
    
  }


  # Boss
  {
    # Run Boss on truth R0
    Ghat_boss_0 <- boss_py(R0, asy.n)
    
    # Run Boss across clouds
    Ghat_boss_unifcloud_001 <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_001[,,i], asy.n))
    Ghat_boss_unifcloud_005 <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_005[,,i], asy.n))
    Ghat_boss_unifcloud_01  <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_01[,,i],  asy.n))
    Ghat_boss_unifcloud_02  <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_02[,,i],  asy.n))
    Ghat_boss_bootstrap     <- lapply(1:ITER, function(i) boss_py(Rhat.bootstrap[,,i],     asy.n))
    
    # CDA metrics — adjacency
    perf_adj.boss_unifcloud_001 <- sapply(Ghat_boss_unifcloud_001, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
    
    perf_adj.boss_unifcloud_005 <- sapply(Ghat_boss_unifcloud_005, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
    
    perf_adj.boss_unifcloud_01 <- sapply(Ghat_boss_unifcloud_01, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
    
    perf_adj.boss_unifcloud_02 <- sapply(Ghat_boss_unifcloud_02, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
    
    perf_adj.boss_bootstrap <- sapply(Ghat_boss_bootstrap, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
    
    # CDA metrics — orientation
    perf_ort.boss_unifcloud_001 <- sapply(Ghat_boss_unifcloud_001, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    
    perf_ort.boss_unifcloud_005 <- sapply(Ghat_boss_unifcloud_005, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    
    perf_ort.boss_unifcloud_01 <- sapply(Ghat_boss_unifcloud_01, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    
    perf_ort.boss_unifcloud_02 <- sapply(Ghat_boss_unifcloud_02, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    
    perf_ort.boss_bootstrap <- sapply(Ghat_boss_bootstrap, function(ghat)
      eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    
    # Plotting — identical layout to PC
    {
      df_adjacency_accuracy_boss <- bind_rows(
        unifcloud_001 = data.frame(value = perf_adj.boss_unifcloud_001[1,], group = "unifcloud_0.01"),
        unifcloud_005 = data.frame(value = perf_adj.boss_unifcloud_005[1,], group = "unifcloud_0.05"),
        unifcloud_01  = data.frame(value = perf_adj.boss_unifcloud_01[1,],  group = "unifcloud_0.1"),
        unifcloud_02  = data.frame(value = perf_adj.boss_unifcloud_02[1,],  group = "unifcloud_0.2"),
        bootstrap     = data.frame(value = perf_adj.boss_bootstrap[1,],     group = "bootstrap")
      )
      
      df_adjacency_sensitivity_boss <- bind_rows(
        unifcloud_001 = data.frame(value = perf_adj.boss_unifcloud_001[2,], group = "unifcloud_0.01"),
        unifcloud_005 = data.frame(value = perf_adj.boss_unifcloud_005[2,], group = "unifcloud_0.05"),
        unifcloud_01  = data.frame(value = perf_adj.boss_unifcloud_01[2,],  group = "unifcloud_0.1"),
        unifcloud_02  = data.frame(value = perf_adj.boss_unifcloud_02[2,],  group = "unifcloud_0.2"),
        bootstrap     = data.frame(value = perf_adj.boss_bootstrap[2,],     group = "bootstrap")
      )
      
      df_orientation_accuracy_boss <- bind_rows(
        unifcloud_001 = data.frame(value = perf_ort.boss_unifcloud_001[1,], group = "unifcloud_0.01"),
        unifcloud_005 = data.frame(value = perf_ort.boss_unifcloud_005[1,], group = "unifcloud_0.05"),
        unifcloud_01  = data.frame(value = perf_ort.boss_unifcloud_01[1,],  group = "unifcloud_0.1"),
        unifcloud_02  = data.frame(value = perf_ort.boss_unifcloud_02[1,],  group = "unifcloud_0.2"),
        bootstrap     = data.frame(value = perf_ort.boss_bootstrap[1,],     group = "bootstrap")
      )
      
      df_orientation_sensitivity_boss <- bind_rows(
        unifcloud_001 = data.frame(value = perf_ort.boss_unifcloud_001[2,], group = "unifcloud_0.01"),
        unifcloud_005 = data.frame(value = perf_ort.boss_unifcloud_005[2,], group = "unifcloud_0.05"),
        unifcloud_01  = data.frame(value = perf_ort.boss_unifcloud_01[2,],  group = "unifcloud_0.1"),
        unifcloud_02  = data.frame(value = perf_ort.boss_unifcloud_02[2,],  group = "unifcloud_0.2"),
        bootstrap     = data.frame(value = perf_ort.boss_bootstrap[2,],     group = "bootstrap")
      )
      
      g_adj_acc_boss <- ggplot(df_adjacency_accuracy_boss, aes(x = group, y = value)) +
        geom_boxplot() + ggtitle("BOSS: Adjacency Accuracy") + theme_minimal() +
        theme(axis.title.x = element_blank())+
        my_scale
      
      g_adj_sen_boss <- ggplot(df_adjacency_sensitivity_boss, aes(x = group, y = value)) +
        geom_boxplot() + ggtitle("BOSS: Adjacency Sensitivity") + theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      
      g_ort_acc_boss <- ggplot(df_orientation_accuracy_boss, aes(x = group, y = value)) +
        geom_boxplot() + ggtitle("BOSS: Orientation Accuracy") + theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      
      g_ort_sen_boss <- ggplot(df_orientation_sensitivity_boss, aes(x = group, y = value)) +
        geom_boxplot() + ggtitle("BOSS: Orientation Sensitivity") + theme_minimal()+
        theme(axis.title.x = element_blank())+
        my_scale
      
      ggboss <- ggpubr::ggarrange(
        g_adj_acc_boss, g_adj_sen_boss,
        g_ort_acc_boss, g_ort_sen_boss,
        ncol = 2, nrow = 2
      )
    }
  }
  
  
}

ggpc
ggboss

Ghat.bootstrap.mean<- Reduce(`+`, Ghat.bootstrap)/ITER
Ghat.unifcloud_005.mean<- Reduce(`+`, Ghat.unifcloud_005)/ITER

round(Ghat.bootstrap.mean,2)
round(Ghat.unifcloud_005.mean,2)
sum(abs(G0- Ghat.bootstrap.mean))
sum(abs(G0- Ghat.unifcloud_005.mean))


if (require(Rgraphviz)) {
  par(mfrow=c(1,4))
  plot(as(t(G0), "graphNEL"), main = "True DAG")
  plot(as(t(Ghat.waldcloud[[1]]), "graphNEL"), main= "Est1")
}  



# Post Check --------------------------------------------------------------
## 1: frob norm of each Rhat with R0
## 2: wald-distance of each Rhat with R0
dwald.waldcloud<- sapply(1:ITER, function(i) wald.test(Rpop= R0, Rsample= Rhat.waldcloud[,,i]))
dwald.unifcloud<- sapply(1:ITER, function(i) wald.test(Rpop= R0, Rsample= Rhat.unifcloud[,,i]))
dwald.corX<- sapply(1:ITER, function(i) wald.test(Rpop= R0, Rsample= Rhat.corX[,,i]))

dwald.waldcloud[5,] %>% mean
dwald.unifcloud[5,] %>% mean
dwald.corX[5,] %>% mean

df <- bind_rows(
  data.frame(value = dwald.waldcloud[1,], group = "waldcloud"),
  data.frame(value = dwald.unifcloud[1,], group = "unifcloud"),
  data.frame(value = dwald.corX[1,], group = "corData")
)

df_chisq<- d*(d-1)/2
ggplot(df, aes(x = value, fill = group)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    bins = 30,
    alpha = 0.35
  ) +
  stat_function(
    fun = dchisq,
    args = list(df = df_chisq),
    color = "black",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  theme_minimal()


# To show that volume ratio is uniformly distributed
ps<- 10
ps2<- ps*(ps-1)/2
ITER= 30000
SAMPLE<- runifcloud(diag(ps),n = asy.n, B=ITER)
SAMPLE[,,1]
dwald.unifcloud<- sapply(1:ITER, function(i) wald.test(Rpop= diag(ps), Rsample= SAMPLE[,,i]))
V =((dwald.unifcloud[1,]/qchisq(0.95, df= ps2))^(ps2/2))
hist(V, freq = F)


# Want to compare clouds
# Truth Layer:
p= 3L; ad= 2L; asy.n= 500L; ITER=2000L
Target = er_dag_py(p=p, ad=ad, n= asy.n, K=ITER)  # Already consists of all simdata set
G0 <- Target$G
R0 <- Target$R
alpha1<- 0.15   # alpha1= 1 minus coveragelevel
alpha2<- 0.05   # used in CI test
if (require(Rgraphviz)) {
  #par(mfrow=c(1,4))
  plot(as(t(G0), "graphNEL"), main = "True DAG")
  #plot(as(t(Ghat.waldcloud[[1]]), "graphNEL"), main= "Est1")
}  


# Returns Array
Rhat.waldcloud<- rwaldcloud(Target$R, n = asy.n, B = ITER)
Rhat.unifcloud<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  alpha1)
Rhat.corX <- lapply(1:ITER, function(i) cor(Target$X[i,,])) %>% simplify2array

# ---- 3D correlation-cloud plot (Plotly) ----
{
  library(plotly)
  
  # ---- helpers ----
  corr_vec <- function(R) {
    c(r12 = R[1, 2], r13 = R[1, 3], r23 = R[2, 3])
  }
  
  cloud_df <- function(A, label) {
    stopifnot(length(dim(A)) == 3L, dim(A)[1] == 3L, dim(A)[2] == 3L)
    B <- dim(A)[3]
    M <- t(vapply(seq_len(B), function(b) corr_vec(A[, , b]), numeric(3)))
    df <- as.data.frame(M)
    df$set <- label
    df$idx <- seq_len(B)
    df
  }
  
  # ---- data ----
  df_corX <- cloud_df(Rhat.corX,      "Bootstrap (corX)")
  df_wald <- cloud_df(Rhat.waldcloud, "Wald cloud")
  df_CI   <- cloud_df(Rhat.unifcloud,   "CI cloud")
  
  R0_vec <- corr_vec(R0)
  df_R0 <- data.frame(
    r12 = R0_vec["r12"],
    r13 = R0_vec["r13"],
    r23 = R0_vec["r23"],
    set = "R0",
    idx = 1
  )
  
  # ---- plot ----
  p <- plot_ly()
  
  # Bootstrap cloud
  p <- add_trace(
    p,
    data = df_corX,
    type = "scatter3d",
    mode = "markers",
    x = ~r12, y = ~r13, z = ~r23,
    name = "R̂ (bootstrap)",
    legendgroup = "bootstrap",
    showlegend = TRUE,
    marker = list(size = 3, opacity = 0.5),
    hoverinfo = "text",
    text = ~paste0(
      "Bootstrap<br>",
      "idx = ", idx,
      "<br>r12 = ", round(r12, 3),
      "<br>r13 = ", round(r13, 3),
      "<br>r23 = ", round(r23, 3)
    )
  )
  
  # Wald cloud
  p <- add_trace(
    p,
    data = df_wald,
    type = "scatter3d",
    mode = "markers",
    x = ~r12, y = ~r13, z = ~r23,
    name = "R̂ (Wald)",
    legendgroup = "wald",
    showlegend = TRUE,
    marker = list(size = 3, opacity = 0.5),
    hoverinfo = "text",
    text = ~paste0(
      "Wald cloud<br>",
      "idx = ", idx,
      "<br>r12 = ", round(r12, 3),
      "<br>r13 = ", round(r13, 3),
      "<br>r23 = ", round(r23, 3)
    )
  )
  
  # CI cloud
  p <- add_trace(
    p,
    data = df_CI,
    type = "scatter3d",
    mode = "markers",
    x = ~r12, y = ~r13, z = ~r23,
    name = "R̂ (CI)",
    legendgroup = "ci",
    showlegend = TRUE,
    marker = list(size = 3, opacity = 0.5),
    hoverinfo = "text",
    text = ~paste0(
      "CI cloud<br>",
      "idx = ", idx,
      "<br>r12 = ", round(r12, 3),
      "<br>r13 = ", round(r13, 3),
      "<br>r23 = ", round(r23, 3)
    )
  )
  
  # True R0
  p <- add_trace(
    p,
    data = df_R0,
    type = "scatter3d",
    mode = "markers",
    x = ~r12, y = ~r13, z = ~r23,
    name = "R₀ (true)",
    legendgroup = "R0",
    showlegend = TRUE,
    marker = list(size = 8, symbol = "diamond", opacity = 1),
    hoverinfo = "text",
    text = ~paste0(
      "R0 (true)<br>",
      "r12 = ", round(r12, 3),
      "<br>r13 = ", round(r13, 3),
      "<br>r23 = ", round(r23, 3)
    )
  )
  
  # ---- layout ----
  p <- layout(
    p,
    title = "3D correlation-space clouds",
    legend = list(
      orientation = "h",
      x = 0,
      y = 1.05,
      title = list(text = "")
    ),
    scene = list(
      xaxis = list(title = "r12 = R[1,2]", range = c(-1, 1)),
      yaxis = list(title = "r13 = R[1,3]", range = c(-1, 1)),
      zaxis = list(title = "r23 = R[2,3]", range = c(-1, 1))
    )
  )
  
  p
  
  
  
  
}

cov.corX= cov(df_corX[,1:3])* asy.n
cov.wald= cov(df_wald[,1:3])* asy.n
#cov(df_CI[,1:3])* asy.n
cov0= metaSEM::asyCov(R0, n= 1)

D= cov.corX- cov0
dd= D[upper.tri(D)]          
sum(dd^2)

D= cov.wald- cov0
dd= D[upper.tri(D)]          
sum(dd^2)


colMeans(df_corX[,1:3])
colMeans(df_wald[,1:3])



# Full Simulation ---------------------------------------------------------

nTarget = 30
p= 12L; ad= 3L; asy.n= 20000L; ITER=300L
ps= p*(p-1)/2
alpha2<- 0.05   # used in CI test

perf_adj.bootstrap<- perf_adj.unifcloud_02 <- perf_adj.unifcloud_01 <- perf_adj.unifcloud_005 <- perf_adj.unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
perf_ort.bootstrap<- perf_ort.unifcloud_02 <- perf_ort.unifcloud_01 <- perf_ort.unifcloud_005 <- perf_ort.unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
perf_adj.boss_bootstrap <- perf_adj.boss_unifcloud_02 <- perf_adj.boss_unifcloud_01 <- perf_adj.boss_unifcloud_005 <- perf_adj.boss_unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
perf_ort.boss_bootstrap <- perf_ort.boss_unifcloud_02 <- perf_ort.boss_unifcloud_01 <- perf_ort.boss_unifcloud_005 <- perf_ort.boss_unifcloud_001 <-array(NA, dim=c(nTarget, 2, ITER))
# For objects above: First row is accuracy, Second row is sensitivity

# SIMULATION Starts here
{
  for(nT in 1:nTarget){
    cat(Sys.time(),";Iteration:", nT, "\n")
    # Truth Layer:
    Target = er_dag_py(p=p, ad=ad, n= asy.n, K=ITER)  # Consists of K data sets
    G0 <- Target$G
    R0 <- Target$R
    X<- Target$X[1,,]   # The only observed data
    
    # Generate Rsamples. Approaches: 1) Bootstrap, 2) 1-alpha uniform wald cloud
    # Returns 3d-Array
    #Rhat.waldcloud<- rwaldcloud(Target$R, n = asy.n, B = ITER)
    #Rhat.corX <- lapply(1:ITER, function(i) cor(Target$X[nT,,])) %>% simplify2array
    Rhat.bootstrap <- lapply(1:ITER, function(i){
      BOOT.ID<- sample(asy.n, asy.n, replace= TRUE)
      X.boot<- X[BOOT.ID,]
      return(cor(X.boot))
    }) %>% simplify2array
    Rhat.unifcloud_001<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.01)
    Rhat.unifcloud_005<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.05)
    Rhat.unifcloud_01<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.1)
    Rhat.unifcloud_02<- runifcloud(Target$R, n = asy.n, B = ITER, alpha =  0.2)
    
    # Run CDA
    # PC
    {
      fit0<- pcalg::pc(suffStat= list(C = R0, n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2 )
      Ghat.0 <- as(fit0, "matrix")
      #fit.waldcloud<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.waldcloud[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      #fit.corX<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.corX[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.unifcloud_001<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_001[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.unifcloud_005<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_005[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.unifcloud_01<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_01[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.unifcloud_02<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.unifcloud_02[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      fit.bootstrap<- lapply(1:ITER, function(i) pcalg::pc(suffStat= list(C = Rhat.bootstrap[,,i], n = asy.n), indepTest = gaussCItest, p= p, alpha = alpha2))
      
      #Ghat.waldcloud<- lapply(fit.waldcloud, function(Y) as(Y, "matrix"))
      #Ghat.corX<- lapply(fit.corX, function(Y) as(Y, "matrix"))
      Ghat.unifcloud_001<- lapply(fit.unifcloud_001, function(Y) as(Y, "matrix"))
      Ghat.unifcloud_005<- lapply(fit.unifcloud_005, function(Y) as(Y, "matrix"))
      Ghat.unifcloud_01<- lapply(fit.unifcloud_01, function(Y) as(Y, "matrix"))
      Ghat.unifcloud_02<- lapply(fit.unifcloud_02, function(Y) as(Y, "matrix"))
      Ghat.bootstrap<- lapply(fit.bootstrap, function(Y) as(Y, "matrix"))
      
      # Compare Ghats with G0
      #perf.waldcloud<- sapply(Ghat.waldcloud, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric)) 
      #perf.corX<- sapply(Ghat.corX, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric)) 
      perf_adj.unifcloud_001[nT,,]<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_adj.unifcloud_005[nT,,]<- sapply(Ghat.unifcloud_005, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_adj.unifcloud_01[nT,,]<- sapply(Ghat.unifcloud_01, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_adj.unifcloud_02[nT,,]<- sapply(Ghat.unifcloud_02, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      perf_adj.bootstrap[nT,,]<- sapply(Ghat.bootstrap, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.adj)) 
      
      perf_ort.unifcloud_001[nT,,]<- sapply(Ghat.unifcloud_001, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.unifcloud_005[nT,,]<- sapply(Ghat.unifcloud_005, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.unifcloud_01[nT,,]<- sapply(Ghat.unifcloud_01, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.unifcloud_02[nT,,]<- sapply(Ghat.unifcloud_02, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
      perf_ort.bootstrap[nT,,]<- sapply(Ghat.bootstrap, function(ghat) (eval_cda(true_adj=t(G0), est_adj=t(ghat))$metric.orient)) 
    }
    
    
    # Boss
    {
      # Run Boss on truth R0
      Ghat_boss_0 <- boss_py(R0, asy.n)
      
      # Run Boss across clouds
      Ghat_boss_unifcloud_001 <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_001[,,i], asy.n))
      Ghat_boss_unifcloud_005 <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_005[,,i], asy.n))
      Ghat_boss_unifcloud_01  <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_01[,,i],  asy.n))
      Ghat_boss_unifcloud_02  <- lapply(1:ITER, function(i) boss_py(Rhat.unifcloud_02[,,i],  asy.n))
      Ghat_boss_bootstrap     <- lapply(1:ITER, function(i) boss_py(Rhat.bootstrap[,,i],     asy.n))
      
      # CDA metrics — adjacency
      perf_adj.boss_unifcloud_001[nT,,] <- sapply(Ghat_boss_unifcloud_001, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      
      perf_adj.boss_unifcloud_005[nT,,] <- sapply(Ghat_boss_unifcloud_005, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      
      perf_adj.boss_unifcloud_01[nT,,] <- sapply(Ghat_boss_unifcloud_01, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      
      perf_adj.boss_unifcloud_02[nT,,] <- sapply(Ghat_boss_unifcloud_02, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      
      perf_adj.boss_bootstrap[nT,,] <- sapply(Ghat_boss_bootstrap, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.adj)
      
      # CDA metrics — orientation
      perf_ort.boss_unifcloud_001[nT,,] <- sapply(Ghat_boss_unifcloud_001, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      
      perf_ort.boss_unifcloud_005[nT,,] <- sapply(Ghat_boss_unifcloud_005, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      
      perf_ort.boss_unifcloud_01[nT,,] <- sapply(Ghat_boss_unifcloud_01, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      
      perf_ort.boss_unifcloud_02[nT,,] <- sapply(Ghat_boss_unifcloud_02, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
      
      perf_ort.boss_bootstrap[nT,,] <- sapply(Ghat_boss_bootstrap, function(ghat)
        eval_cda(true_adj = t(G0), est_adj = t(ghat))$metric.orient)
    }
  }
}

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

perf_ort.pc
perf_adj.boss
perf_ort.boss
