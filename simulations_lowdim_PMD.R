packages <- c("foreach", "doSNOW", "RobSparseMVA", "PMA")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
ncl <- 100

simulation_function <- function(r, scenario, setting, contamination) {
  set.seed(r * 10)
  sample <- Sample(scenario = scenario)
  data <- simulate_Sample(sample, setting = setting, cont_ratio = contamination)
  true_CCA <- get_true_cca(sample)
  rank <- 2
  summary <- c()

  try(
    {
      t0 <- Sys.time()
      CCA_PMA_perm <- PMA::CCA.permute(data$x, data$y)
      res_CCA <- PMA::CCA(data$x, data$y,
                          penaltyx = CCA_PMA_perm$bestpenaltyx,
                          penaltyz = CCA_PMA_perm$bestpenaltyz,
                          K = rank)

      t1 <- Sys.time()

      angle_a <- RobSparseMVA::principal_angles(true_CCA$a, res_CCA$u)$angles
      angle_b <- RobSparseMVA::principal_angles(true_CCA$b, res_CCA$v)$angles
      TPR_a <- RobSparseMVA::TPR(as.matrix(true_CCA$a[, 1:rank]), as.matrix(res_CCA$u))
      TPR_b <- RobSparseMVA::TPR(as.matrix(true_CCA$b[, 1:rank]), as.matrix(res_CCA$v))
      TNR_a <- RobSparseMVA::TNR(true_CCA$a[, 1:rank], res_CCA$u)
      TNR_b <- RobSparseMVA::TNR(true_CCA$b[, 1:rank], res_CCA$v)

      summary <- list(
        angle_a = angle_a, angle_b = angle_b,
        a = res_CCA$u, b = res_CCA$v,
        TPR_a = TPR_a, TPR_b = TPR_b,
        TNR_a = TNR_a, TNR_b = TNR_b,
        pen_x = res_CCA$penaltyx,
        pen_y = res_CCA$penaltyz,
        corr = res_CCA$cors,
        runtime = t1 - t0
      )
    },
    FALSE
  )

  return(summary)
}


R <- 100
simestcl <- makeCluster(ncl)
clusterCall(simestcl, function() { source("simulation_design.R") })
registerDoSNOW(simestcl)

Low_HO <- foreach(
  r = 1:R,
  .packages = c("PMA","RobSparseMVA")
) %dopar% {
  res_normal <- simulation_function(r, "LOW_HO", "NORMAL", 0)
  res_cont <- simulation_function(r, "LOW_HO", "CONTAMINATED", 0.05)
  res_tdist <- simulation_function(r, "LOW_HO", "T_DIST", 0)

  list(
    res_normal = res_normal,
    res_cont = res_cont,
    res_tdist = res_tdist
  )
}
save(Low_HO, file = file.path("results", "low_ho_PMD.Rdata"))

stopCluster(simestcl)
