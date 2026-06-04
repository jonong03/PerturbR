# script to run in msi

library(future.apply)
library(data.table)

source("docs/WorkingCode/Functions.R")   # source required functions
source("UseCase/msi/param_setup.R")      # source parameters for simulation

args <- commandArgs(trailingOnly = TRUE)
j <- as.integer(args[1])

dir.create("UseCase/msi/sim_outputs", showWarnings = FALSE)

fout <- sprintf("UseCase/msi/sim_outputs/res_param_%03d.rds", j)

if (file.exists(fout)) {
  message("Output already exists for j = ", j, ". Skipping.")
  quit(save = "no")
}

source("docs/WorkingCode/Functions.R")
message("Running parameter set ", j, "/", nrow(param_grid))

rho_j  <- rho_[param_grid$rho_id[j]]
beta_j <- beta_[[param_grid$beta_id[j]]]
ss_j   <- ss_[param_grid$ss_id[j]]

M_j <- gen_AR1(rho = rho_j, p = nvar)

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

plan(multisession, workers = n_workers)
runs_j <- future_lapply(seq_len(B), function(b) {
  library(mvtnorm)
  simdriver(Rxx = M_j, beta = beta_j, k = ktrue, ss = ss_j, nchain = nchain)
}, future.seed = TRUE)

out.point_j <- array(NA_real_, dim = c(1, B))
out.chain_j <- array(NA_real_, dim = c(nchain, B))
out.boot_j  <- array(NA_real_, dim = c(nchain, B))

for (i in seq_len(B)) {
  out.point_j[1, i] <- runs_j[[i]]$khat.point
  out.chain_j[, i] <- runs_j[[i]]$khat.chain
  out.boot_j[, i]  <- runs_j[[i]]$khat.boot
}

saveRDS(
  list(
    j = j,
    param = param_grid[j, ],
    out.point = out.point_j,
    out.chain = out.chain_j,
    out.boot  = out.boot_j
  ),
  file = fout
)

plan(sequential)

message("Saved ", fout)