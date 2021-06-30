library(vegan)
library(ggplot2)


# function to get PCO scores in data frame, # of positive and negative eigenvalues
# input should be a matrix of distance/dissimilarity
# and a data frame of factors (factor_unit) with "year" and "site" columns
# input_matrix and factor_df should have the same number of rows (samples)

get_pco_scores <- function(input_matrix, factor_df) {
  
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


# function to get dt and dbt for multivariate control charts
# input should be the pco_data from get_pco_scores() 
# and b (the number of time points to be used as baseline for dbt)

mcc <- function(pco_data, b) {
  
  pco_df <- pco_data$pco_df
  p <- pco_data$p
  n <- pco_data$n
  
  site <- unique(pco_df$site)
  year <- unique(pco_df$year)
  
  dt_output <- data.frame(year = character(0), site = character(0), dt = numeric(0))
  dbt_output <- data.frame(year = character(0), site = character(0), dbt = numeric(0))
  
  for (i in 1:length(site)) {
    
    subset_site <- site[i]
    dat <- pco_df[pco_df$site == subset_site, ]
    rownames(dat) <- dat$year
    # get squared Euclidean dist. among years
    dat_sq <- as.matrix((dist(dat[, 1:p])) ^ 2 - (dist(dat[, (p+1):(p+n)])) ^ 2)
    
    cc_dat_sum <- data.frame(year = colnames(dat_sq), t = 1:ncol(dat_sq), 
                             sst = NA, dt = NA, 
                             ssbt = NA, dbt = NA)
    
    # dt
    for (j in 1:nrow(cc_dat_sum)) {
      # the sum of squared interpoint dissimilarities among all points up to time t divided by t
      cc_dat_sum$sst[j] <- (sum(dat_sq[1:j, 1:j]) / 2) / cc_dat_sum$t[j]
    }
    for (k in 2:nrow(cc_dat_sum)) {
      cc_dat_sum$dt[k] <- sqrt(cc_dat_sum$t[k] / (cc_dat_sum$t[k-1]) * (cc_dat_sum$sst[k] - cc_dat_sum$sst[k-1]))
    }
    
    # dbt
    ssb <- cc_dat_sum$sst[cc_dat_sum$t == b]
    for (l in (b+1):nrow(cc_dat_sum)) {
      cc_dat_sum$ssbt[l] <- ((sum(dat_sq[1:b, 1:b]) / 2) + sum(dat_sq[l, 1:b])) / (b + 1)
      cc_dat_sum$dbt[l] <- sqrt((b + 1) / b * (cc_dat_sum$ssbt[l] - ssb))
    }
    
    dt_temp <- data.frame(year = cc_dat_sum$year, site = rep(subset_site, nrow(cc_dat_sum)), dt = cc_dat_sum$dt)
    dt_output <- rbind(dt_output, dt_temp)
    
    dbt_temp <- data.frame(year = cc_dat_sum$year, site = rep(subset_site, nrow(cc_dat_sum)), dbt = cc_dat_sum$dbt)
    dbt_output <- rbind(dbt_output, dbt_temp)
    
  }
  
  mcc_data <- list(dt_output, dbt_output, site, year)
  names(mcc_data) <- c("dt_output", "dbt_output", 'site', 'year')
  
  return(mcc_data)
  
}


# function for bootstrapping to get confidence intervals
# input should be the pco_data and mcc_data from the above functions
# and b and nboot (the number of bootstrap)

mcc_bootstrap <- function(pco_data, mcc_data, b, nboot) {
  
  pco_df <- pco_data$pco_df
  p <- pco_data$p
  n <- pco_data$n
  
  dt_output <- mcc_data$dt_output
  dbt_output <- mcc_data$dbt_output
  site <- mcc_data$site
  year <- mcc_data$year
  
  dt_boot <- matrix(nrow = nrow(dt_output), ncol = nboot)
  dbt_boot <- matrix(nrow = nrow(dbt_output), ncol = nboot)
  
  for(i in 1:nboot) {
    
    dt_output_boot <- data.frame(year = character(0), site = character(0), dt = numeric(0))
    dbt_output_boot <- data.frame(year = character(0), site = character(0), dbt = numeric(0))
    
    for (j in 1:length(site)) {
      
      subset_site <- site[j]
      subdat <- pco_df[pco_df$site == subset_site, ]
      
      vec <- 1:nrow(subdat)
      boot_vec <- sample(vec, replace = TRUE)
      dat_boot <- subdat[boot_vec, 1:(p+n)]
      rownames(dat_boot) <- subdat$year
      
      dat_sq_boot <- as.matrix((dist(dat_boot[, 1:p])) ^ 2 - (dist(dat_boot[, (p+1):(p+n)])) ^ 2)
      
      cc_dat_sum_boot <- data.frame(year = colnames(dat_sq_boot), t = 1:ncol(dat_sq_boot), 
                                    sst = NA, dt = NA, 
                                    ssbt = NA, dbt = NA)
      
      # dt boot
      for (k in 1:nrow(cc_dat_sum_boot)) {
        cc_dat_sum_boot$sst[k] <- (sum(dat_sq_boot[1:k, 1:k]) / 2) / cc_dat_sum_boot$t[k]
      }
      for (l in 2:nrow(cc_dat_sum_boot)) {
        cc_dat_sum_boot$dt[l] <- sqrt(cc_dat_sum_boot$t[l] / (cc_dat_sum_boot$t[l-1]) * (cc_dat_sum_boot$sst[l] - cc_dat_sum_boot$sst[l-1]))
      }
      
      # dbt boot
      ssb_boot <- cc_dat_sum_boot$sst[cc_dat_sum_boot$t == b]
      for (m in (b+1):nrow(cc_dat_sum_boot)) {
        cc_dat_sum_boot$ssbt[m] <- ((sum(dat_sq_boot[1:b, 1:b]) / 2) + sum(dat_sq_boot[m, 1:b])) / (b + 1)
        cc_dat_sum_boot$dbt[m] <- sqrt((b + 1) / b * (cc_dat_sum_boot$ssbt[m] - ssb_boot))
      }
      
      dt_temp_boot <- data.frame(year = cc_dat_sum_boot$year, site = rep(subset_site, nrow(cc_dat_sum_boot)), dt = cc_dat_sum_boot$dt)
      dt_output_boot <- rbind(dt_output_boot, dt_temp_boot)
      
      dbt_temp_boot <- data.frame(year = cc_dat_sum_boot$year, site = rep(subset_site, nrow(cc_dat_sum_boot)), dbt = cc_dat_sum_boot$dbt)
      dbt_output_boot <- rbind(dbt_output_boot, dbt_temp_boot)
      
    }
    
    dt_boot[, i] <- dt_output_boot$dt
    dbt_boot[, i] <- dbt_output_boot$dbt
  }
  
  dt_boot_df <- as.data.frame(dt_boot)
  dt_data <- cbind(dt_output, dt_boot_df)
  
  dbt_boot_df <- as.data.frame(dbt_boot)
  dbt_data <- cbind(dbt_output, dbt_boot_df)
  
  mcc_boot_data <- list(dt_data, dbt_data)
  names(mcc_boot_data) <- c("dt_data", "dbt_data")
  
  return(mcc_boot_data)
  
}


### MCC and bootstrap, continued from 01_data_prep.R

# pick either fish or coral data

load("rda_files/kala_fish.rda")
# load("rda_files/kaho_coral_mcc.rda")

# get Bray-Curtis dissimilarity (fish or coral)

# fish
fish_sqrt <- sqrt(fish_unit_mcc)   # sqrt transform fish count
fish_bray <- as.matrix(vegdist(fish_sqrt, "bray")) * 100  # Bray-Curtis dissimilarity scaled from 0 to 100

# coral
# coral_bray <- as.matrix(vegdist(coral_unit_mcc, "bray")) * 100

pco_data <- get_pco_scores(fish_bray, factor_unit)  # input should be fish_bray or coral_bray
save(pco_data, file = "rda_files/kala_fish_pco.rda")  # name appropriately

mcc_data <- mcc(pco_data, b = 2)

mcc_boot_data <- mcc_bootstrap(pco_data, mcc_data, b = 2, nboot = 1000)

save(mcc_boot_data, file = "rda_files/kala_fish_mcc.rda")  # give an appropriate name

rm(list = ls())

