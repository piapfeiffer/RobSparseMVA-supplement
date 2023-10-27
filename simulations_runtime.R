packages <- c("foreach", "doSNOW", "RobSparseMVA")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
ncl <- 10

simulation_function <- function(r, scenario, setting, methods, contamination, kwargs) {
  set.seed(r * 10)
  sample <- Sample(scenario = scenario)
  data <- simulate_Sample(sample,
    setting = setting,
    cont_ratio = contamination,
    cont_strength = 2
  )
  true_CCA <- get_true_cca(sample)
  rank <- 1
  res_summary <- list()

  for (method in methods) {
    try({
    data_x <- (data$x - matrix(apply(data$x, 2, median), byrow = TRUE, ncol = ncol(data$x), nrow = nrow(data$x)))
    data_y <- (data$y - matrix(apply(data$y, 2, median), byrow = TRUE, ncol = ncol(data$y), nrow = nrow(data$y)))

    t0 <- Sys.time()
    res_CCA <- RobSparseMVA::ccaMM(data_x, data_y,
      method = method,
      kwargs[[method]],
      nearPD = FALSE,
      alpha_x = rep(1, rank),
      alpha_y = rep(1, rank),
      k = rank,
      tol = 1e-3, lr = 1e-2,
      epochs = 300,
      criterion = "TPO",
      penalties = list(
        pen_x = 1.104315,
        pen_y = 1.104315
      )
    )
    t1 <- Sys.time()

    angle_a <- RobSparseMVA::principal_angles(true_CCA$a, res_CCA$a)$angles
    angle_b <- RobSparseMVA::principal_angles(true_CCA$b, res_CCA$b)$angles
    TPR_a <- RobSparseMVA::TPR(as.matrix(true_CCA$a[, 1:rank]), as.matrix(res_CCA$a))
    TPR_b <- RobSparseMVA::TPR(as.matrix(true_CCA$b[, 1:rank]), as.matrix(res_CCA$b))
    TNR_a <- RobSparseMVA::TNR(true_CCA$a[, 1:rank], res_CCA$a)
    TNR_b <- RobSparseMVA::TNR(true_CCA$b[, 1:rank], res_CCA$b)

    res_summary[[method]] <- list(
      angle_a = angle_a, angle_b = angle_b,
      a = res_CCA$a, b = res_CCA$b,
      TPR_a = TPR_a, TPR_b = TPR_b,
      TNR_a = TNR_a, TNR_b = TNR_b,
      pen_x = res_CCA$pen_x,
      pen_y = res_CCA$pen_y,
      corr = res_CCA$measure,
      runtime = t1 - t0
    )
    }, FALSE)
  }

  return(res_summary)
}


R <- 2 # 10
simestcl <- makeCluster(ncl)
clusterCall(simestcl, function() {
  source("simulation_design.R")
})
registerDoSNOW(simestcl)

Runtime_Dim <- foreach(
  r = 1:R,
  .packages = c(
    "RobSparseMVA"
  ),
  .errorhandling = "pass"
) %dopar% {
  methods <- c("Pearson", "Spearman", "Kendall", "MRCD", "OGK")
  kwargs <- list(
    Pearson = NA,
    Spearman = NA,
    Kendall = NA,
    MRCD = c(alpha = 0.75),
    OGK = NA
  )
  res_50_100 <- simulation_function(r, "RUNTIME_50_100", "NORMAL", methods, 0, kwargs)
  res_100_100 <- simulation_function(r, "RUNTIME_100_100", "NORMAL", methods, 0, kwargs)
  res_500_100 <- simulation_function(r, "RUNTIME_500_100", "NORMAL", methods, 0, kwargs)
  res_1000_100 <- simulation_function(r, "RUNTIME_1000_100", "NORMAL", methods, 0, kwargs)
  res_5000_100 <- simulation_function(r, "RUNTIME_5000_100", "NORMAL", methods, 0, kwargs)
  res_10000_100 <- simulation_function(r, "RUNTIME_10000_100", "NORMAL", methods, 0, kwargs)

  list(
    res_50_100 = res_50_100,
    res_100_100 = res_100_100,
    res_500_100 = res_500_100,
    res_1000_100 = res_1000_100,
    res_5000_100 = res_5000_100,
    res_10000_100 = res_10000_100
  )
}
save(Runtime_Dim, file = file.path("results", "runtime_dim_rank.Rdata"))

stopCluster(simestcl)
