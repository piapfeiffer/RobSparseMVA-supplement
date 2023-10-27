# File to replicate results for precision evaluations
# ===================================================

library(RobSparseMVA)

# high-dimensional setting
set.seed(10)
sample <- Sample(scenario = "HIGH_HO")
data <- simulate(sample, setting = "NORMAL")
true_CCA <- get_true_cca(sample)

p <- ncol(data$x)
q <- ncol(data$y)
n <- nrow(data$x)

k <- 2

CCA_obj <- CCA("Pearson", alpha_x = rep(1, k), alpha_y = rep(1, k))

A <- matrix(NA, nrow = p, ncol = k)
B <- matrix(NA, nrow = q, ncol = k)

PHI <- matrix(NA, nrow = n, ncol = k)
ETA <- matrix(NA, nrow = n, ncol = k)

CORR <- rep(NA, k)
PEN_X <- rep(NA, k)
PEN_Y <- rep(NA, k)

for (i in 1:k) {
  PEN_X[i] <- norm(as.matrix(true_CCA$a[,i]), "1")
  PEN_Y[i] <- norm(as.matrix(true_CCA$b[,i]), "1")
}

SUMMARY <- list()


C <- rbind(cbind(sample$cov_xx, sample$cov_xy), cbind(Matrix::t(sample$cov_xy), sample$cov_yy))

# true covariance, theoretically optimal penalties
for (i in 1:k){
  res <- train_CCA(C, p, q,
                   PEN_X[i], PEN_Y[i] , i,
                   CCA_obj$alpha_x[i], CCA_obj$alpha_y[i],
                   1, tol = 1e-5, lr = 1e-2, 1000,
                   low_a = as.matrix(A[,1:max(1,(i-1))]),
                   low_b = as.matrix(B[,1:max(1,(i-1))]),
                   use_transform = FALSE,
                   init_ortho = TRUE)

  A[, i] <- res$best_a
  B[, i] <- res$best_b

  PHI[, i] <- data_x %*% res$best_a
  ETA[, i] <- data_y %*% res$best_b

  CORR[i] <- res$best_measure
}

angle_a <- principal_angles(true_CCA$a, A)$angles
angle_b <- principal_angles(true_CCA$b, B)$angles

TPR_a <- TPR(true_CCA$a[,2], A[,2])
TPR_b <- TPR(true_CCA$b[,2], B[,2])
TNR_a <- TNR(true_CCA$a[,2], A[,2])
TNR_b <- TNR(true_CCA$b[,2], B[,2])

# compute hyperparameters via BO
A_hp <- matrix(NA, nrow = ncol(data_x), ncol = k)
B_hp <- matrix(NA, nrow = ncol(data_y), ncol = k)

CORR_hp <- rep(NA, k)
PEN_X_hp <- rep(NA, k)
PEN_Y_hp <- rep(NA, k)

for (i in 1:k) {
    res_param <- bayesian_optimization_CCA(C, p, q, n, CCA_obj, i,
                                           low_a = as.matrix(true_CCA$a[,1:max(1,(i-1))]),
                                           low_b = as.matrix(true_CCA$b[,1:max(1,(i-1))]),
                                           tol = 1e-5, lr = 1e-2, 1000,
                                           criterion = "TPO",
                                           init_ortho = TRUE)

    PEN_X_hp[i] <- res_param$best_params$pen_x
    PEN_Y_hp[i] <- res_param$best_params$pen_y

    SUMMARY[[i]] <- res_param$summary

  res <- train_CCA(C, p, q,
                   PEN_X_hp[i] , PEN_Y_hp[i] , i,
                   CCA_obj$alpha_x[i], CCA_obj$alpha_y[i],
                   1, tol = 1e-5, lr = 1e-2, 1000,
                   low_a = as.matrix(A_hp[,1:max(1,(i-1))]),
                   low_b = as.matrix(B_hp[,1:max(1,(i-1))]),
                   init_ortho = TRUE)

  A_hp[, i] <- res$best_a
  B_hp[, i] <- res$best_b

  CORR_hp[i] <- res$best_measure

}

angle_a_hp <- principal_angles(true_CCA$a, A_hp)$angles
angle_b_hp <- principal_angles(true_CCA$b, B_hp)$angles

TPR_a_hp <- TPR(true_CCA$a[,2], A_hp[,2])
TPR_b_hp <- TPR(true_CCA$b[,2], B_hp[,2])
TNR_a_hp <- TNR(true_CCA$a[,2], A_hp[,2])
TNR_b_hp <- TNR(true_CCA$b[,2], B_hp[,2])

# high-dimensional setting
set.seed(10)
sample <- Sample(scenario = "HIGH_HO")
data <- simulate(sample, setting = "NORMAL")
true_CCA <- get_true_cca(sample)

p <- ncol(data$x)
q <- ncol(data$y)
n <- nrow(data$x)

k <- 2

CCA_obj <- CCA("Pearson", alpha_x = rep(1, k), alpha_y = rep(1, k))

A <- matrix(NA, nrow = p, ncol = k)
B <- matrix(NA, nrow = q, ncol = k)

PHI <- matrix(NA, nrow = n, ncol = k)
ETA <- matrix(NA, nrow = n, ncol = k)

CORR <- rep(NA, k)
PEN_X <- rep(NA, k)
PEN_Y <- rep(NA, k)

for (i in 1:k) {
  PEN_X[i] <- norm(as.matrix(true_CCA$a[,i]), "1")
  PEN_Y[i] <- norm(as.matrix(true_CCA$b[,i]), "1")
}

SUMMARY <- list()


C <- rbind(cbind(sample$cov_xx, sample$cov_xy), cbind(Matrix::t(sample$cov_xy), sample$cov_yy))

# true covariance, theoretically optimal penalties
for (i in 1:k){
  res <- train_CCA(C, p, q,
                   PEN_X[i], PEN_Y[i] , i,
                   CCA_obj$alpha_x[i], CCA_obj$alpha_y[i],
                   1, tol = 1e-5, lr = 1e-2, 1000,
                   low_a = as.matrix(A[,1:max(1,(i-1))]),
                   low_b = as.matrix(B[,1:max(1,(i-1))]),
                   use_transform = FALSE,
                   init_ortho = TRUE)

  A[, i] <- res$best_a
  B[, i] <- res$best_b

  PHI[, i] <- data_x %*% res$best_a
  ETA[, i] <- data_y %*% res$best_b

  CORR[i] <- res$best_measure
}

angle_a <- principal_angles(true_CCA$a, A)$angles
angle_b <- principal_angles(true_CCA$b, B)$angles

TPR_a <- TPR(true_CCA$a[,2], A[,2])
TPR_b <- TPR(true_CCA$b[,2], B[,2])
TNR_a <- TNR(true_CCA$a[,2], A[,2])
TNR_b <- TNR(true_CCA$b[,2], B[,2])

# compute hyperparameters via BO
A_hp <- matrix(NA, nrow = ncol(data_x), ncol = k)
B_hp <- matrix(NA, nrow = ncol(data_y), ncol = k)

CORR_hp <- rep(NA, k)
PEN_X_hp <- rep(NA, k)
PEN_Y_hp <- rep(NA, k)

for (i in 1:k) {
    res_param <- bayesian_optimization_CCA(C, p, q, n, CCA_obj, i,
                                           low_a = as.matrix(true_CCA$a[,1:max(1,(i-1))]),
                                           low_b = as.matrix(true_CCA$b[,1:max(1,(i-1))]),
                                           tol = 1e-5, lr = 1e-2, 1000,
                                           criterion = "TPO",
                                           init_ortho = TRUE)

    PEN_X_hp[i] <- res_param$best_params$pen_x
    PEN_Y_hp[i] <- res_param$best_params$pen_y

    SUMMARY[[i]] <- res_param$summary

  res <- train_CCA(C, p, q,
                   PEN_X_hp[i] , PEN_Y_hp[i] , i,
                   CCA_obj$alpha_x[i], CCA_obj$alpha_y[i],
                   1, tol = 1e-5, lr = 1e-2, 1000,
                   low_a = as.matrix(A_hp[,1:max(1,(i-1))]),
                   low_b = as.matrix(B_hp[,1:max(1,(i-1))]),
                   init_ortho = TRUE)

  A_hp[, i] <- res$best_a
  B_hp[, i] <- res$best_b

  CORR_hp[i] <- res$best_measure

}

angle_a_hp <- principal_angles(true_CCA$a, A_hp)$angles
angle_b_hp <- principal_angles(true_CCA$b, B_hp)$angles

TPR_a_hp <- TPR(true_CCA$a[,2], A_hp[,2])
TPR_b_hp <- TPR(true_CCA$b[,2], B_hp[,2])
TNR_a_hp <- TNR(true_CCA$a[,2], A_hp[,2])
TNR_b_hp <- TNR(true_CCA$b[,2], B_hp[,2])