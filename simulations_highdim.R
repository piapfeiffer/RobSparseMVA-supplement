# File to replicate results for high-dimensional setting
# ======================================================
packages <- c("foreach", "doSNOW", "RobSparseMVA")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
ncl <- 100

simulation_function <- function(r,
                                scenario,
                                setting,
                                methods,
                                contamination,
                                kwargs) {
  set.seed(r * 10)
  sample <- Sample(scenario = scenario)
  data <- simulate_Sample(sample,
    setting = setting,
    cont_ratio = contamination,
    cont_strength = 2
  )
  true_CCA <- get_true_cca(sample)
  order <- 2
  res_summary <- list()

  for (method in methods) {
    try(
      {
        data_x <- (data$x - matrix(apply(data$x, 2, median),
          byrow = TRUE,
          ncol = ncol(data$x),
          nrow = nrow(data$x)
        ))
        data_y <- (data$y - matrix(apply(data$y, 2, median),
          byrow = TRUE,
          ncol = ncol(data$y),
          nrow = nrow(data$y)
        ))

        t0 <- Sys.time()
        res_CCA <- RobSparseMVA::ccaMM(data_x, data_y,
          method = method,
          kwargs[[method]],
          nearPD = TRUE,
          alpha_x = rep(1, order),
          alpha_y = rep(1, order),
          k = order,
          tol = 1e-5, lr = 1e-3,
          epochs = 300,
          criterion = "TPO"
        )


        t1 <- Sys.time()

        angle_a <- RobSparseMVA::principal_angles(true_CCA$a, res_CCA$a)$angles
        angle_b <- RobSparseMVA::principal_angles(true_CCA$b, res_CCA$b)$angles
        TPR_a <- RobSparseMVA::TPR(as.matrix(true_CCA$a[, 1:order]), 
                                   as.matrix(res_CCA$a))
        TPR_b <- RobSparseMVA::TPR(as.matrix(true_CCA$b[, 1:order]), 
                                   as.matrix(res_CCA$b))
        TNR_a <- RobSparseMVA::TNR(true_CCA$a[, 1:order], res_CCA$a)
        TNR_b <- RobSparseMVA::TNR(true_CCA$b[, 1:order], res_CCA$b)

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
      },
      TRUE
    )
  }

  return(res_summary)
}


R <- 100
simestcl <- makeCluster(ncl)
pb <- txtProgressBar(min = 1, max = R, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
clusterCall(simestcl, function() {
  source("simulation_design.R")
})
registerDoSNOW(simestcl)

High_HO <- foreach(
  r = 1:R,
  .packages = c("RobSparseMVA"),
  .errorhandling = "pass",
  .options.snow = opts
) %dopar% {
  methods <- c("Pearson", "Spearman", "Kendall", "MRCD", "OGK")
  kwargs <- list(Pearson = NA, 
                 Spearman = NA, 
                 Kendall = NA, 
                 MRCD = c(alpha = 0.75), 
                 OGK = NA)
  res_normal <- simulation_function(r, "HIGH_HO", "NORMAL", methods, 0, kwargs)
  res_cont <- simulation_function(r, "HIGH_HO", "CONTAMINATED", methods, 0.05, kwargs)
  res_tdist <- simulation_function(r, "HIGH_HO", "T_DIST", methods, 0, kwargs)

  setTxtProgressBar(pb, r)

  list(
    res_normal = res_normal,
    res_cont = res_cont,
    res_tdist = res_tdist
  )
}
save(High_HO, file = file.path("results", "high_ho_orthstart.Rdata"))

stopCluster(simestcl)
