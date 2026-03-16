# Inspecting BOSS outcome

rm(list=ls()); gc()
pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, reticulate, pcalg, ggplot2, ggpubr)
use_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/.venv/bin/python", required = TRUE)
py_config() 
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/process1.py")
reticulate::source_python("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR-CDA/boss_py.py")
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/PerturbR/docs/WorkingCode/Functions.R", local = TRUE)
source("/Users/jonong/Library/CloudStorage/OneDrive-Personal/Documents/1- Projects/CDA/Scripts/Functions.R", local = TRUE)


## Predefined functions
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

# 1: Generate trueR
set.seed(123)
p=5L; ad=2L; asy.n= 50L; B= 500
df<- m<- p*(p-1)/2  # degree of freedom
Target <- er_dag_py(p = p, ad = ad, n = asy.n, K = 500L)
G0 <- Target$G
R0 <- Target$R

# 3: Learn Boss

R1 = cor(Target$X[89,,])   #id=2, 89
pc = run_pc_mat(R1, asy.n, alpha = 0.05)
boss = boss_py(R1, asy.n)

fit =  pcalg::pc(
  suffStat = list(C = R1, n = asy.n),
  indepTest = gaussCItest,
  p = p,
  alpha = 0.05
)
as(fit, "amat")

pc
boss

get_metric(G0 = G0, Ghat= pc)
get_metric(G0 = G0, Ghat= boss)

par(mfrow=c(1,3))
if(require(Rgraphviz)){
  plot(as(G0, "graphNEL"),main="Truth")
  plot(fit, main="PC")
  plot(as(boss, "graphNEL"), main="BOSS")
} 
eval_cda(true_adj= G0, pc)
eval_cda(G0, boss)



# 2: Make perturbations
Rhat<- rwaldcloud(R0= R0, n= asy.n, B=B)
#pcl<- lapply(1:B, function(i) run_pc_mat(Rhat[,,i], n_obs= asy.n, alpha = 0.05))
#bossl <- lapply(1:B, function(i) boss_py(Rhat[,,i], asy.n))

pcl<- lapply(1:B, function(i) run_pc_mat(cor(Target$X[i,,]), n_obs= asy.n, alpha = 0.05))
bossl <- lapply(1:B, function(i) boss_py(cor(Target$X[i,,]), asy.n))

outpc<- lapply(pcl, function(Ghat) get_metric(G0, Ghat))
outpc.adj<- sapply(outpc, function(Q) Q$adj) %>% t()
outpc.ort<- sapply(outpc, function(Q) Q$orient) %>% t()

outboss<- lapply(bossl, function(Ghat) get_metric(G0, Ghat))
outboss.adj<- sapply(outboss, function(Q) Q$adj) %>% t()
outboss.ort<- sapply(outboss, function(Q) Q$orient) %>% t()

par(mfrow=c(1,2))
hist(outpc.adj)
hist(outboss.adj)
summary(outpc.adj)
summary(outboss.adj)
which.max(outpc.adj[,1]- outboss.adj[,1])

?pc
