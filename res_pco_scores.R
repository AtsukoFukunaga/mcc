# Function to get PCO scores after covariates are fitted
# Input should be a matrix of distance/dissimilarity
# and a data frame of factors (factor_df).
# Include "latitude" and "longitude" columns in the factor_df for the location variables 
# to be fitted as a PCA score.
# Also specify the names of columns to be used as covariate.
# input_matrix and factor_df should have the same number of rows (samples).

get_res_pco_scores <- function(input_matrix, factor_df, factor_names) {
  
  require(vegan)
  
  # calculate matrix A
  a_bray <- -0.5 * input_matrix * input_matrix  
  # get identity matrix
  In <- diag(nrow(a_bray))
  # get matrix of 1/#samples
  ln <- matrix(1 / nrow(a_bray), nrow = nrow(a_bray), ncol = ncol(a_bray))
  # get Gower's centered matrix
  g_bray <- (In - ln) %*% a_bray %*% (In - ln)
  # eigenvalue decomposition
  g_eig <- eigen(g_bray)
  # scale eigenvectors by sqrt(eigenvalues) and get principal coordinates
  g_pco_score <- g_eig$vectors %*% diag(abs(g_eig$values) ^ 0.5)
  
  # covariate
  if ("latitude" %in% colnames(factor_df) & "longitude" %in% colnames(factor_df)) {
    loc <- prcomp(factor_df[, c("latitude","longitude")], center = T, scale = F)
    pc1 <- predict(loc, newdat = factor_df[, c("latitude","longitude")])[, 1]
    pc1_sc <- scale(pc1)
    
    cov <- as.data.frame(pc1_sc)
    cov <- cbind(cov, factor_df[, factor_names])
    colnames(cov) <- c("loc", factor_names)
  } else {
    cov <- factor_df[, factor_names]
    colnames(cov) <- factor_names
  }
  xcr <- model.matrix(~.-1, data = cov)
  
  # fitted
  h <- xcr %*% (solve(t(xcr) %*% xcr)) %*% t(xcr)
  outpro_fit <- h %*% g_bray %*% h
  dbrda <- eigen(outpro_fit)
  dbscore <- dbrda$vectors %*% diag(abs(dbrda$values) ^ 0.5)
  r_square <- sum(diag(outpro_fit))/sum(diag(g_bray))
  cat("R-square value for the fitted model: ", r_square)

  # residuals
  i <- diag(nrow(g_bray))
  res <- (i - h) %*% g_bray %*% (i - h)
  res_eig <- eigen(res)
  res_pco_score <- res_eig$vectors %*% diag(abs(res_eig$values) ^ 0.5)
  
  p <- length(res_eig$values[res_eig$values >= 0]) ## # of positive eigenvalues
  n <- length(res_eig$values[res_eig$values < 0]) ## # of negative eigenvalues
  
  pco_df <- cbind(as.data.frame(res_pco_score), factor_df)
  
  pco_data <- list(pco_df, p, n)
  names(pco_data) <- c("pco_df", "p", "n")
  
  return(pco_data)
}
