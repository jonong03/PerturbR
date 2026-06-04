# Sim 5: Reverse Simulation

rm(list=ls());gc()
source("docs/WorkingCode/Functions.R")

pacman::p_load(metaSEM, dplyr, data.table, mvtnorm, rethinking, future.apply, parallel, parallelly, ggplot2)

d=30; n= 1e3; B=1000
Rhat<- rlkjcorr(1,K=d,eta=0.1)
kappa(metaSEM::asyCov(Rhat, n= n))
Rs<- rwaldcloud(Rhat, n= n, B=B)   # Sample 1000 R* around Rhat
out<- sapply(1:B, function(b) wald.test(Rs[,,b], Rhat, fisherz= FALSE, alpha = 0.05)) %>% t() %>% as.data.table(.)   # Check if R*'s CI contains Rhat. If reject=1, no containment
mean(out$pval < 0.05, na.rm=TRUE)


# SIM FUNCTION
sim_reverse_containment<- function(d, n, B, eta){
  
  is_psd <- function(A, tol = 1e-4) {
    A <- as.matrix(A)
    min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) >= tol
  }
  
  psd_check= FALSE
  i=0; max_iter= 100
  while(psd_check!= TRUE){
    i= i+1
    Rhat<- rlkjcorr(1,K=d,eta=eta)
    PSI<- metaSEM::asyCov(Rhat, n=n)
    psd_check<- is_psd(PSI)
    
    if(i> max_iter){break}
  }
  
  Rs<- rwaldcloud(Rhat, n= n, B=B)   # Sample 1000 R* around Rhat
  
  # Check if R*'s CI contains Rhat. 
  # alpha here does not matter, interest is to get p-value
  out<- sapply(1:B, function(b) wald.test(Rs[,,b], Rhat, fisherz= FALSE, alpha = 0.05)) %>%
    t() %>% 
    as.data.table(.)   
  
  return(out$pval)
  
}

sim_reverse_containment(d=3, n=100000, B=1000, eta=5)

d_list<- c(5, 10,20,30)
n_list <- c(100, 200, 500,
            1e3, 2e3, 5e3,
            1e4, 2e4, 5e4,
            1e5, 2e5, 5e5,
            1e6)
B<- 300
eta_list <- c(10, 5, 1, 0.5)

param_grid<- as.data.table(expand.grid(d= d_list, n= n_list, B= B, eta= eta_list))
param_grid


# parallel run

plan(multisession)
pvals_list <- future_lapply(seq_len(nrow(param_grid)), function(i) {
  tryCatch(
    sim_reverse_containment(param_grid$d[i], param_grid$n[i], param_grid$B[i], param_grid$eta[i]),
    error = function(e) e
  )
}, future.seed = TRUE)


cutoff= 0.1
param_grid[, `:=`(
  pvals = lapply(pvals_list, function(x) if (inherits(x, "error")) NA_real_ else x),
  err  = lapply(pvals_list, function(x) if (inherits(x, "error")) conditionMessage(x) else NA_character_)
)][, `:=`(
  contain_rate = sapply(pvals, function(x) mean(x > cutoff, na.rm = TRUE)),
  mean_pval   = sapply(pvals, function(x) mean(x, na.rm = TRUE)),
  na_prop     = sapply(pvals, function(x) mean(is.na(x))),
  pval_med = sapply(pvals, median, na.rm = TRUE),
  pval_q25 = sapply(pvals, quantile, probs = 0.25, na.rm = TRUE),
  pval_q75 = sapply(pvals, quantile, probs = 0.75, na.rm = TRUE)
)]

g_contain<- ggplot(param_grid,
                    aes(x = n, y = contain_rate, color = factor(eta), group = eta)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~ d, labeller = labeller(d = function(x) paste0("#variables = ", x))) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x="Sample size (n, log10)", y="Reverse containment rate", color=expression(eta),
       title= "Reverse Containment Rate",
       subtitle= paste("containment: when p-val >", cutoff)) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

g_medianp<- ggplot(param_grid, aes(n, pval_med, color = factor(eta), fill = factor(eta), group = eta)) +
  geom_ribbon(aes(ymin = pval_q25, ymax = pval_q75), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ d, labeller = labeller(d = function(x) paste0("#variables = ", x))) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x="Sample size (n, log10)", y="Median p-value (IQR ribbon)", color=expression(eta), fill=expression(eta),
       title= "Median (with IQR) of p-value",
       subtitle= ""
       ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

g_naprop<- ggplot(param_grid, aes(x = as.factor((n)), y = factor(eta), fill = na_prop)) +
  geom_tile(color = "grey80", linewidth = 0.3) +
  facet_wrap(~ d, labeller = labeller(d = function(x) paste0("#variables = ", x))) +
  labs(x = "log10(n)", y = expression(eta), fill = "NA proportion",
       title= "Wald Cloud Breakdown Map",
       subtitle= "High NA proportion implies the wald geometry region is ill-defined"
       ) +
  scale_fill_viridis_c(option = "A")+
  ggtitle("Wald Cloud Breakdown Map")+
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )


g_contain
g_medianp
g_naprop
# Asymptotic variance is unstable, not PSD, R lies near boundary, a map region where the wald geometry is ill defined



B=1000
Rb<- rwaldcloud(R1, n= 1e5, B=B, output="full")
out<- sapply(1:B, function(b) wald.test(R1, Rb[,,b], asy.n=1e5)) %>% t()
mean(out[,5])

# Check theoretical and actual density

R_d3<- rlkjcorr(1,K=3, eta=1)
R_d5<- rlkjcorr(1, K=5, eta=1)
R_d10<- rlkjcorr(1,K=10, eta=1)
R_d20<- rlkjcorr(1, K=20, eta=1)

## Simulate actual and theoretical density of test statistic
density.simulate<- function(R0, B=1000, asy.n=100000){
  d= ncol(R0)
  p=d*(d-1)/2
  Rb<- rwaldcloud(R0, n= asy.n, B=B, output="full")
  out<- sapply(1:B, function(b) wald.test(Rpop= R0, Rsample= Rb[,,b], asy.n= asy.n))
  x<- out[1,] # test statistic
  
  hist(x, breaks = 30, probability = TRUE, main = paste0("#Nodes (d):",d, "; degree-of-freedom (p)=",p) ,xlab = "Wald-Statistic")
  curve(dchisq(x, df= p),  add = TRUE, col = "red", lwd = 2)
}

par(mfrow=c(2,2), oma = c(0, 0, 3, 0))

density.simulate(R_d3)
density.simulate(R_d5)
density.simulate(R_d10)
density.simulate(R_d20)

mtext("Empirical histogram with Chi-square curve (red)", outer = TRUE, cex = 1.4, line = 1)


