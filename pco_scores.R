# Function to get PCO scores in data frame, # of positive and negative eigenvalues
# Input should be a matrix of distance/dissimilarity
# and a data frame of factors (factor_unit) with "year" and "site" columns.
# input_matrix and factor_df should have the same number of rows (samples).

get_pco_scores <- function(input_matrix, factor_df) {
  
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
  
  pco_df <- cbind(as.data.frame(g_pco_score), factor_df)
  
  p <- length(g_eig$values[g_eig$values >= 0])  # of positive eigenvalues
  n <- length(g_eig$values[g_eig$values < 0])  # of negative eigenvalues
  
  pco_data <- list(pco_df, p, n)
  names(pco_data) <- c("pco_df", "p", "n")
  
  return(pco_data)
  
}
