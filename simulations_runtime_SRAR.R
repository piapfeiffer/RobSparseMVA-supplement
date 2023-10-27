packages <- c("foreach", "doSNOW", "RobSparseMVA")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
ncl <- 10

simulation_function <- function(r, scenario, setting, contamination) {
  set.seed(r * 10)
  sample <- Sample(scenario = scenario)
  data <- simulate_Sample(sample, setting = setting, cont_ratio = contamination)
  true_CCA <- get_true_cca(sample)
  rank <- 1
  summary <- c()

  try(
    {
      data_x <- (data$x-matrix(apply(data$x,2,median),byrow=TRUE,ncol=ncol(data$x),nrow=nrow(data$x)))
      data_y <- (data$y-matrix(apply(data$y,2,median),byrow=TRUE,ncol=ncol(data$y),nrow=nrow(data$y)))
      t0 <- Sys.time()
      res_CCA <- RobSparseMVA::ccaAR(data_x, data_y, rank = rank,
                                               lambdaAseq=0.1104315 / ncol(data_x),
                                               lambdaBseq=0.1104315 / ncol(data_y))
      t1 <- Sys.time()

      angle_a <- RobSparseMVA::principal_angles(true_CCA$a, res_CCA$ALPHA)$angles
      angle_b <- RobSparseMVA::principal_angles(true_CCA$b, res_CCA$BETA)$angles
      TPR_a <- RobSparseMVA::TPR(as.matrix(true_CCA$a[, 1:rank]), as.matrix(res_CCA$ALPHA))
      TPR_b <- RobSparseMVA::TPR(as.matrix(true_CCA$b[, 1:rank]), as.matrix(res_CCA$BETA))
      TNR_a <- RobSparseMVA::TNR(true_CCA$a[, 1:rank], res_CCA$ALPHA)
      TNR_b <- RobSparseMVA::TNR(true_CCA$b[, 1:rank], res_CCA$BETA)

      summary <- list(
        angle_a = angle_a, angle_b = angle_b,
        a = res_CCA$ALPHA, b = res_CCA$BETA,
        TPR_a = TPR_a, TPR_b = TPR_b,
        TNR_a = TNR_a, TNR_b = TNR_b,
        pen_x = res_CCA$lambdaA,
        pen_y = res_CCA$lambdaB,
        corr = res_CCA$cancors,
        runtime = t1 - t0
      )
  },
  FALSE
  )

  return(summary)
}


R <- 10
simestcl <- makeCluster(ncl)
clusterCall(simestcl, function() {
  source("simulation_design.R")
})
registerDoSNOW(simestcl)

Runtime_Dim_SRAR <- foreach(
  r = 1:R,
  .packages = c(
    "RobSparseMVA"
  )
) %dopar% {
  res_50_100 <- simulation_function(r, "RUNTIME_50_100", "NORMAL", 0)
  res_100_100 <- simulation_function(r, "RUNTIME_100_100", "NORMAL", 0)
  res_500_100 <- simulation_function(r, "RUNTIME_500_100", "NORMAL", 0)
  res_1000_100 <- simulation_function(r, "RUNTIME_1000_100", "NORMAL", 0)
  res_5000_100 <- simulation_function(r, "RUNTIME_5000_100", "NORMAL", 0)
  res_10000_100 <- simulation_function(r, "RUNTIME_10000_100", "NORMAL", 0)

  list(res_50_100 = res_50_100,
       res_100_100 = res_100_100,
       res_500_100 = res_500_100,
       res_1000_100 = res_1000_100,
       res_5000_100 = res_5000_100,
       res_10000_100 = res_10000_100)
}
save(Runtime_Dim_SRAR, file = file.path("results", "runtime_dim_SRAR.Rdata"))

stopCluster(simestcl)
