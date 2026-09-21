## Precompute the (p0, p1) grid of variance estimators used by the
## "Differences in the Variance estimator" section of the
## Binomial_and_Count_Outcomes vignette.
##
## Run this script to regenerate inst/vignettes/BinomOutcomes_vals.rds.

library(SteppedPower)
library(parallel)

icc    <- 0.2
p_vals <- seq(0.02, 0.98, by = 0.01)
DM     <- construct_DesMat(Cl = rep(1,4), N = 50)
p_grid <- expand.grid(p0 = p_vals, p1 = p_vals)

## Flatten the (p0, p1) grid into independent cells
one_cell <- function(p0, p1) {
  p_mid <- (p0 + p1) / 2
  tau_re <- sqrt(p_mid * (1 - p_mid) * icc / (1 - icc))

  re <- suppressMessages(glsPower(DesMat = DM, mu0 = p0, mu1 = p1,
               tau = tau_re, family = "binomial",
               verbose = 2)$VarianceMatrix[1, 1])
  al <- suppressMessages(glsPower(DesMat = DM, mu0 = p0, mu1 = p1,
               alpha_0_1_2 = c(icc, icc),
               family = "binomial",
               verbose = 2)$VarianceMatrix[1, 1])
  lin <- suppressMessages(glsPower(DesMat = DM, mu0 = p0, mu1 = p1,
               tau = tau_re, sigma = sqrt(p_mid * (1 - p_mid)),
               family = "gaussian",
               verbose = 2)$VarianceMatrix[1, 1])
  c(re = re, al = al, lin = lin)
}

ncores <- max(1, parallel::detectCores() - 1, na.rm = TRUE)
cl     <- parallel::makeCluster(ncores)

parallel::clusterExport(cl, c("one_cell", "p_grid", "icc","DM"), envir = environment())
parallel::clusterEvalQ(cl, library(SteppedPower)) |> invisible()

res <- parallel::parLapply(cl, seq_len(nrow(p_grid)), function(k) {
  one_cell(p_grid$p0[k], p_grid$p1[k])
})
parallel::stopCluster(cl)

vals <- do.call(rbind, res)
out  <- cbind(p_grid,vals)


out_path <- file.path("inst", "vignettes", "BinomOutcomes_vals.rds")
saveRDS(out, out_path)
