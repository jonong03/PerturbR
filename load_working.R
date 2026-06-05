files <- list.files("UseCase/msi/sim_outputs", pattern = "\\.rds$", full.names = TRUE)
length(files)

out.point <- array(NA_real_, dim = c(1, B, length(files)))
out.boot  <- array(NA_real_, dim = c(nchain, B, length(files)))
out.chain <- array(NA_real_, dim = c(nchain, B, length(files)))

for (f in files) {
  x <- readRDS(f)
  j <- x$j
  
  out.point[1, , j] <- x$out.point[1,]
  out.chain[, , j] <- x$out.chain
  out.boot[, , j]  <- x$out.boot
}


hist(out.point)


CI.mcmc<- apply(out.chain, 2, range)
CI.boot<- apply(out.boot, 2, function(x) quantile(x, c(0.025, 0.975)))
ktrue=0.9

mean((CI.mcmc[1,] < ktrue) & (CI.mcmc[2,] > ktrue))
mean((CI.boot[1,] < ktrue) & (CI.boot[2,] > ktrue))

