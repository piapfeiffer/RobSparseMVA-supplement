load("data/data_FTIR_HOG.Rdata")

packages <- c("foreach", "doSNOW", "RobSparseMVA")
suppressMessages(suppressWarnings(packages <- lapply(packages, FUN = function(x) {
  library(x, character.only = TRUE)
})))
ncl <- 5

R <- 5
simestcl <- makeCluster(ncl)
registerDoSNOW(simestcl)
parallel::clusterExport(simestcl, "data")

pb <- txtProgressBar(min = 1, max = R, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

CV_residuals_P <- foreach(
  r = 1:R,
  .packages = c(
    "RobSparseMVA"
  ),
  .errorhandling = "pass",
  .options.snow = opts
) %dopar% {
  set.seed(10*r)
  
  n_train <- base::sample(1:nrow(data$FTIR), floor(0.9 * nrow(data$FTIR)), replace = FALSE)
  n_test <- c(1:nrow(data$FTIR))[-n_train]
  
  rem_FTIR <- which(apply(data$FTIR[n_train,], 2, mad) == 0)
  FTIR <- data$FTIR[,-rem_FTIR]
  rem_HOG <- which(apply(data$HOG[n_train,], 2, mad) == 0)
  HOG <- data$HOG[,-rem_HOG]
  
  med_HOG <- apply(HOG[n_train,], 2, median)
  mad_HOG <- apply(HOG[n_train,], 2, mad)
  centered_HOG <- HOG[n_train,] - med_HOG
  scaled_HOG <- (HOG[n_train,] - med_HOG) / mad_HOG
  
  med_FTIR <- apply(FTIR[n_train,], 2, median)
  mad_FTIR <- apply(FTIR[n_train,], 2, mad)
  centered_FTIR <- FTIR[n_train,] - med_FTIR
  scaled_FTIR <- (FTIR[n_train,] - med_FTIR) / mad_FTIR
  
  res_P <- RobSparseMVA::ccaMM(FTIR[n_train,], scaled_HOG, 
                               method = "Pearson",
                               k = 1, 
                               nearPD = FALSE,
                               alpha_x = rep(1, 1),
                               alpha_y = rep(0, 1),
                               tol = 1e-3, lr = 1e-2, epochs = 2000,
                               lr_decay = 0.99,
                               criterion = "TPO")
  CV_P <-  ((FTIR[n_test,]) %*% res_P$a - sign(res_P$measure) * ((HOG[n_test,]-med_HOG)/mad_HOG) %*% res_P$b)^2
  setTxtProgressBar(pb, r) 
  list(n_train = n_train,
       residual_P = CV_P,
       res_P = res_P)
}
save(CV_residuals_P, file = file.path("results_tribology","cv_residuals_P_FTIR_HOGscale.Rdata"))

CV_residuals_OGK <- foreach(
  r = 1:R,
  .packages = c(
    "RobSparseMVA",
  ),
  .errorhandling = "pass",
  .options.snow = opts
) %dopar% {
  n_train <- base::sample(1:nrow(data$FTIR), floor(0.9 * nrow(data$FTIR)), replace = FALSE)
  n_test <- c(1:nrow(data$FTIR))[-n_train]
  
  rem_FTIR <- which(apply(data$FTIR[n_train,], 2, mad) == 0)
  FTIR <- data$FTIR[,-rem_FTIR]
  rem_HOG <- which(apply(data$HOG[n_train,], 2, mad) == 0)
  HOG <- data$HOG[,-rem_HOG]
  
  med_HOG <- apply(HOG[n_train,], 2, median)
  mad_HOG <- apply(HOG[n_train,], 2, mad)
  centered_HOG <- HOG[n_train,] - med_HOG
  scaled_HOG <- (HOG[n_train,] - med_HOG) / mad_HOG
  
  med_FTIR <- apply(FTIR[n_train,], 2, median)
  mad_FTIR <- apply(FTIR[n_train,], 2, mad)
  centered_FTIR <- FTIR[n_train,] - med_FTIR
  scaled_FTIR <- (FTIR[n_train,] - med_FTIR) / mad_FTIR
  
  res_OGK <- RobSparseMVA::ccaMM(FTIR[n_train,], scaled_HOG, 
                               method = "OGK",
                               k = 1, 
                               nearPD = FALSE,
                               alpha_x = rep(1, 1),
                               alpha_y = rep(0, 1),
                               tol = 1e-3, lr = 1e-2, epochs = 2000,
                               lr_decay = 0.99,
                               criterion = "TPO")
  CV_OGK <-  ((FTIR[n_test,]) %*% res_OGK$a - sign(res_OGK$measure) * ((HOG[n_test,]-med_HOG)/mad_HOG) %*% res_OGK$b)^2
  setTxtProgressBar(pb, r) 
  list(n_train = n_train,
       residual_OGK = CV_OGK,
       res_OGK = res_OGK)
}
save(CV_residuals_OGK, file = file.path("results","cv_residuals_OGK_FTIR_HOGscale.Rdata"))
stopCluster(simestcl)
