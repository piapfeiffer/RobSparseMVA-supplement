packages <- c("foreach", "doSNOW", "RobSparseMVA", "glmnet")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
source(file.path("scca_code", "SCCA-code.R"))

ncl <- 4

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
      SCCA_cv <- cv.SCCA.equal(data$x, data$y,
                               lambda = seq(from = 0.2, to = 0.01, length = 10),
                               npairs = rank,
                               init.method = "svd")
      res_CCA <- SCCA(data$x, data$y,
                      lambda.alpha = SCCA_cv$bestlambda,
                      lambda.beta = SCCA_cv$bestlambda,
                      npairs = rank,
                      init.method = "svd")

      t1 <- Sys.time()

      angle_a <- RobSparseMVA::principal_angles(true_CCA$a, res_CCA$alpha)$angles
      angle_b <- RobSparseMVA::principal_angles(true_CCA$b, res_CCA$beta)$angles
      TPR_a <- RobSparseMVA::TPR(as.matrix(true_CCA$a[, 1:rank]), as.matrix(res_CCA$alpha))
      TPR_b <- RobSparseMVA::TPR(as.matrix(true_CCA$b[, 1:rank]), as.matrix(res_CCA$beta))
      TNR_a <- RobSparseMVA::TNR(true_CCA$a[, 1:rank], res_CCA$alpha)
      TNR_b <- RobSparseMVA::TNR(true_CCA$b[, 1:rank], res_CCA$beta)

      summary <- list(
        angle_a = angle_a, angle_b = angle_b,
        a = res_CCA$alpha, b = res_CCA$beta,
        TPR_a = TPR_a, TPR_b = TPR_b,
        TNR_a = TNR_a, TNR_b = TNR_b,
        pen_x = SCCA_cv$bestlambda,
        pen_y = SCCA_cv$bestlambda,
        corr = diag(cor(data$x %*% res_CCA$alpha, data$y %*% res_CCA$beta)),
        runtime = t1 - t0
      )
    },
    FALSE
  )

  return(summary)
}


R <- 100
simestcl <- makeCluster(ncl)
clusterCall(simestcl, function() { source("simulation_design.R"); source(file.path("scca_code", "SCCA-code.R")) })
registerDoSNOW(simestcl)

Low_HO <- foreach(
  r = 1:R,
  .packages = c("RobSparseMVA", "glmnet")
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

save(Low_HO, file = file.path("results", "low_ho_SCCA.Rdata"))

stopCluster(simestcl)
