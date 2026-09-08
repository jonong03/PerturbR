library(data.table)

B <- 500          # Monte Carlo replicates
BOOT.ITER <- 300  # Bootstrap samples per Monte Carlo replicate
nchain <- 3000    # PerturbR samples per Monte Carlo replicate

R2true <- c(0.5)
phi <- c(0)
ss <- c(300, 500, 1000, 3000, 5000)
p <- c(11)

beta.list <- c(
  #"function(p) rep(0.6, p)",
  "function(p) rep(0.2, p)"
  #"function(p) rep(c(-0.4, 0.4), length = p)"
)

param.m <- expand.grid(
  R2true = R2true,
  phi = phi,
  ss = ss,
  p = p,
  beta = beta.list,
  stringsAsFactors = FALSE
)

param.m <-   as.data.table(param.m)

param.m$beta_value <- Map(
  function(f, p) eval(parse(text = f))(p),
  param.m$beta,
  param.m$p
)

param.m[, scenarioID := .I]


