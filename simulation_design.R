# Simulation Settings
# contains all scenarios that are presented in the associated paper.

#' sample internal constructor
new_Sample <- function(scenario) {
  # constructor for sample class
  if (scenario == "LOW_HO") {
    p <- 10
    q <- 10
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    diag(cov_xx) <- 1
    diag(cov_yy) <- 1
    cov_xy <- matrix(0, nrow = p, ncol = q)
    cov_xy[1, 1] <- 0.9
    cov_xy[2, 2] <- 0.7
  } else if (scenario == "HIGH_HO") {
    p <- 100
    q <- 100
    n <- 50

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.9
    cov_xx[11:20, 11:20] <- 0.7
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.9
    cov_yy[11:20, 11:20] <- 0.7
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.9
    cov_xy[11:20, 11:20] <- 0.5
  } else if (scenario == "RUNTIME_50_100") {
    p <- 10
    q <- 50
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.8
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.8
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.8
  } else if (scenario == "RUNTIME_100_100") {
    p <- 10
    q <- 100
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.8
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.8
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.8
  } else if (scenario == "RUNTIME_500_100") {
    p <- 10
    q <- 500
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.8
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.8
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.8
  } else if (scenario == "RUNTIME_1000_100") {
    p <- 10
    q <- 1000
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.8
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.8
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.8
  } else if (scenario == "RUNTIME_5000_100") {
    p <- 10
    q <- 5000
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.8
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.8
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.8
  } else if (scenario == "RUNTIME_10000_100") {
    p <- 10
    q <- 10000
    n <- 100

    cov_xx <- matrix(0, ncol = p, nrow = p)
    cov_yy <- matrix(0, ncol = q, nrow = q)
    cov_xy <- matrix(0, nrow = p, ncol = q)

    cov_xx[1:10, 1:10] <- 0.8
    diag(cov_xx) <- 1

    cov_yy[1:10, 1:10] <- 0.8
    diag(cov_yy) <- 1

    cov_xy[1:10, 1:10] <- 0.8
  }

  structure(
    list(
      scenario = scenario,
      n = n,
      p = p,
      q = q,
      cov_xx = cov_xx,
      cov_yy = cov_yy,
      cov_xy = cov_xy
    ),
    class = "Sample"
  )
}

#' constructor for Sample
#' @export
Sample <- function(scenario) {
  # user-friendly helper for sample class
  new_Sample(scenario)
}

#' Simulate sample from Sample class
#' @export
simulate_Sample <- function(object,
                            setting,
                            cont_ratio = 0.1,
                            cont_strength = 2,
                            type = "two-sided") {
  clean <- 1 - cont_ratio
  sigma <- rbind(
    cbind(object$cov_xx, object$cov_xy),
    cbind(Matrix::t(object$cov_xy), object$cov_yy)
  )

  if (setting == "NORMAL") {
    data <- mvtnorm::rmvnorm(floor(object$n),
      mean = rep(0, object$p + object$q),
      sigma = sigma, checkSymmetry = F
    )
  } else if (setting == "CONTAMINATED") {
    cov_xy_cont <- matrix(0, nrow = object$p, ncol = object$q)
    sigma_cont <- rbind(
      cbind(object$cov_xx, cov_xy_cont),
      cbind(t(cov_xy_cont), object$cov_yy)
    )

    if (type == "two-sided") {
      data_clean <- mvtnorm::rmvnorm(floor(clean * object$n),
        mean = rep(0, object$p + object$q),
        sigma = sigma
      )
      data_cont <- mvtnorm::rmvnorm(object$n - floor(clean * object$n),
        mean = rep(cont_strength, object$p + object$q),
        sigma = sigma_cont
      )
    } else if (type == "one-sided") {
      data_clean <- mvtnorm::rmvnorm(floor(object$n),
        mean = rep(0, object$p + object$q),
        sigma = sigma, checkSymmetry = F
      )
      data_cont <- mvtnorm::rmvnorm(object$n - floor(clean * object$n),
        mean = rep(cont_strength, object$p + object$q),
        sigma = sigma_cont
      )

      data_cont <- cbind(data_clean[
        floor(clean * object$n + 1):(floor(object$n)),
        1:object$p
      ], data_cont[, (object$p + 1):(object$p + object$q)])
      data_clean <- data_clean[1:floor(clean * object$n), ]
    }

    data <- rbind(data_clean, data_cont)
  } else if (setting == "T_DIST") {
    df <- 3
    data <- mvtnorm::rmvnorm(floor(object$n),
      mean = rep(0, object$p + object$q),
      sigma = sigma
    ) / matrix(
      rep(
        sqrt(rchisq(object$n,
          df = df,
          ncp = 0
        ) / df),
        object$p + object$q
      ),
      byrow = F,
      ncol = object$p + object$q
    )
  }

  cat("Sample simulated from setting", setting, "\n")
  return(list(
    x = as.matrix(data[, 1:object$p]),
    y = as.matrix(data[, (object$p + 1):(object$p + object$q)])
  ))
}


#' get true CCA sample
#' @export
get_true_cca <- function(object, ...) {
  UseMethod("get_true_cca")
}

#' get true CCA from sample from Sample class
#' @export
get_true_cca.Sample <- function(object, rank, ...) {
  SigmaXX_decomp <- eigen(object$cov_xx)
  SigmaXX_decomp$values <- ifelse(SigmaXX_decomp$values < 0, 0, SigmaXX_decomp$values)
  SigmaXX_negsqrt <- SigmaXX_decomp$vectors %*% diag(1 / sqrt(SigmaXX_decomp$values), length(SigmaXX_decomp$values)) %*% t(SigmaXX_decomp$vectors)
  SigmaYY_decomp <- eigen(object$cov_yy)
  SigmaYY_decomp$values <- ifelse(SigmaYY_decomp$values < 0, 0, SigmaYY_decomp$values)
  SigmaYY_negsqrt <- SigmaYY_decomp$vectors %*% diag(1 / sqrt(SigmaYY_decomp$values), length(SigmaYY_decomp$values)) %*% t(SigmaYY_decomp$vectors)

  trueCCA <- svd(SigmaXX_negsqrt %*% object$cov_xy %*% SigmaYY_negsqrt)
  true_rho <- round(trueCCA$d, 5)
  rank <- length(which(true_rho != 0))
  true_a <- SigmaXX_negsqrt %*% trueCCA$u[, 1:rank]
  true_b <- SigmaYY_negsqrt %*% trueCCA$v[, 1:rank]

  for (i in 1:rank) {
    true_a[, i] <- true_a[, i]
    true_b[, i] <- true_b[, i]
  }
  return(list(
    rho = true_rho,
    rank = rank,
    a = true_a,
    b = true_b
  ))
}
