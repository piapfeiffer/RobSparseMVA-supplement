packages <- c("foreach", "doSNOW", "RobSparseMVA")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
ncl <- 100

simulation_function <- function(r, scenario, setting, methods, contamination, kwargs) {
    set.seed(r * 10)
    sample <- Sample(scenario = scenario)
    data <-simulate_Sample(sample, setting = setting, cont_ratio = contamination)
    true_CCA <- get_true_cca(sample)
    rank <- 2
    summary <- c()

    for (method in methods){

      try({
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
                                     k = rank, 
                                     alpha_x = rep(1, rank),
                                     alpha_y = rep(1, rank),
                                     tol = 1e-3, 
                                     lr = 1e-1, 
                                     epochs = 200)
      t1 <- Sys.time()

      angle_a <- RobSparseMVA::principal_angles(true_CCA$a, res_CCA$a)$angles
      angle_b <- RobSparseMVA::principal_angles(true_CCA$b, res_CCA$b)$angles
      TPR_a <- RobSparseMVA::TPR(as.matrix(true_CCA$a[, 1:rank]), as.matrix(res_CCA$a))
      TPR_b <- RobSparseMVA::TPR(as.matrix(true_CCA$b[, 1:rank]), as.matrix(res_CCA$b))
      TNR_a <- RobSparseMVA::TNR(true_CCA$a[, 1:rank], res_CCA$a)
      TNR_b <- RobSparseMVA::TNR(true_CCA$b[, 1:rank], res_CCA$b)

      summary[[method]] <- list(
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
    return(summary)
}


R <- 100
simestcl <- makeCluster(ncl)
clusterCall(simestcl, function() {
  source("simulation_design.R")
})
registerDoSNOW(simestcl)

Low_HO_break <- foreach(
  r = 1:R,
  .packages = c(
    "RobSparseMVA"
  )
) %dopar% {
  methods <- c("Pearson", "Spearman", "Kendall", "MRCD", "OGK")
  kwargs <- list(Pearson = NA, 
                 Spearman = NA, 
                 Kendall = NA, 
                 MRCD = c(alpha = 0.75), 
                 OGK = NA)
  res_normal <- simulation_function(r, "LOW_HO", "NORMAL", methods, 0, kwargs)
  res_cont1 <- simulation_function(r, "LOW_HO", "CONTAMINATED", methods, 0.1, kwargs)
  res_cont2 <- simulation_function(r, "LOW_HO", "CONTAMINATED", methods, 0.2, kwargs)
  res_cont3 <- simulation_function(r, "LOW_HO", "CONTAMINATED", methods, 0.3, kwargs)
  res_cont4 <- simulation_function(r, "LOW_HO", "CONTAMINATED", methods, 0.4, kwargs)
  res_cont5 <- simulation_function(r, "LOW_HO", "CONTAMINATED", methods, 0.5, kwargs)

  list(res_normal = res_normal,
       res_cont1 = res_cont1,
       res_cont2 = res_cont2,
       res_cont3 = res_cont3,
       res_cont4 = res_cont4,
       res_cont5 = res_cont5)
}
save(Low_HO_break, file = file.path("results","low_ho_break.Rdata"))

stopCluster(simestcl)



