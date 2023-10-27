library(CCA)
library(RobSparseMVA)
library(PMA)

data("nutrimouse")

rank <- 1
rem_gene <- which(apply(nutrimouse$gene, 2, mad) == 0)
if (length(rem_gene) > 0) {
  gene <- as.matrix(nutrimouse$gene[,-rem_gene])
} else {gene <- nutrimouse$gene}
rem_lipid <- which(apply(nutrimouse$lipid, 2, mad) == 0)
if (length(rem_lipid) > 0) {
  lipid <- as.matrix(nutrimouse$lipid[,-rem_lipid])
} else {lipid <- nutrimouse$lipid}

med_gene <- apply(gene, 2, median)
mad_gene <- apply(gene, 2, mad)
centered_gene <- gene - med_gene
scaled_gene <- (gene - med_gene) / mad_gene

med_lipid <- apply(lipid, 2, median)
mad_lipid <- apply(lipid, 2, mad)
centered_lipid <- lipid - med_lipid
scaled_lipid <- (lipid - med_lipid) / mad_lipid

res_OGK <- RobSparseMVA::ccaMM(scaled_gene, scaled_lipid, 
                                   method = "OGK",
                                   k = rank, 
                                   alpha_x = rep(1, rank), 
                                   alpha_y = rep(1, rank),
                                   tol = 1e-5, lr = 1e-3, epochs = 2000,
                                   lr_decay = 1,
                                   criterion = "TPO")

save(res_OGK, file = file.path("results_nutrimouse","res_OGK.Rdata"))


# LOOCV to assess predictive properties
CV_OGK <- rep(NA, nrow(nutrimouse$gene))
CV_PMA <- rep(NA, nrow(nutrimouse$gene))
measure_OGK <- rep(NA, nrow(nutrimouse$gene))
measure_PMA <- rep(NA, nrow(nutrimouse$gene))

for (r in 1:nrow(nutrimouse$gene)){
  set.seed(10*r)
  n_train <- c(1:nrow(nutrimouse$gene))[-r]
  n_test <- c(1:nrow(nutrimouse$gene))[-n_train]
  
  rank <- 1
  rem_gene <- which(apply(nutrimouse$gene[n_train,], 2, mad) == 0)
  if (length(rem_gene) > 0) {
    gene <- as.matrix(nutrimouse$gene[,-rem_gene])
  } else {gene <- nutrimouse$gene}
  rem_lipid <- which(apply(nutrimouse$lipid[n_train,], 2, mad) == 0)
  if (length(rem_lipid) > 0) {
    lipid <- as.matrix(nutrimouse$lipid[,-rem_lipid])
  } else {lipid <- nutrimouse$lipid}
  
  med_gene <- apply(gene[n_train,], 2, median)
  mad_gene <- apply(gene[n_train,], 2, mad)
  centered_gene <- gene[n_train,] - med_gene
  scaled_gene <- (gene[n_train,] - med_gene) / mad_gene
  
  med_lipid <- apply(lipid[n_train,], 2, median)
  mad_lipid <- apply(lipid[n_train,], 2, mad)
  centered_lipid <- lipid[n_train,] - med_lipid
  scaled_lipid <- (lipid[n_train,] - med_lipid) / mad_lipid

    CCA_obj_OGK <- RobSparseMVA::CCA("OGK", alpha_x = rep(1, rank), alpha_y = rep(1, rank))
  res_OGK <- RobSparseMVA::compute_CCA(CCA_obj_OGK, scaled_gene, scaled_lipid, k = rank,
                                     tol = 1e-5, lr = 1e-3, epochs = 2000,
                                     lr_decay = 1,
                                     criterion = "TPO")
  CV_OGK[r] <-  (as.matrix((gene[n_test,]-med_gene)/mad_gene) %*% res_OGK$a - sign(res_OGK$measure) *t(as.matrix((lipid[n_test,]-med_lipid)/mad_lipid)) %*% res_OGK$b)^2
  measure_OGK[r] <- res_OGK$measure
  
  CCA_PMA_perm <- PMA::CCA.permute(scaled_gene, scaled_lipid)
  res_CCA <- PMA::CCA(scaled_gene, scaled_lipid,
                      penaltyx = CCA_PMA_perm$bestpenaltyx,
                      penaltyz = CCA_PMA_perm$bestpenaltyz,
                      K = rank)
  CV_PMA[r] <-  (as.matrix((gene[n_test,]-med_gene)/mad_gene) %*% res_CCA$u - sign(res_CCA$cors) *t(as.matrix((lipid[n_test,]-med_lipid)/mad_lipid)) %*% res_CCA$v)^2
  measure_PMA[r] <- res_CCA$cors
}
save(CV_PMA, file = file.path("results_nutrimouse","cv_residuals_PMA_notsquared.Rdata"))
save(CV_OGK, file = file.path("results_nutrimouse","cv_residuals_OGK_notsquared.Rdata"))
save(measure_PMA, file = file.path("results_nutrimouse","measure_PMA.Rdata"))
save(measure_OGK, file = file.path("results_nutrimouse","measure_OGK.Rdata"))