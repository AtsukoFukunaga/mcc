# Function to get dt and dbt for multivariate control charts.
# Input should be the pco_data from get_pco_scores()
# and b (the number of time points to be used as baseline for dbt).

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


# Function for bootstrapping to get confidence bounds later.
# Input should be the pco_data and mcc_data from the above functions (get_pco_scores() and mcc())
# and b and nboot (the number of bootstrap).

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


# Functions to calculate 50 and 95 percentiles

qt95 <- function(dat) {
  quantile(dat, probs = 0.95, na.rm = T)
}
qt50 <- function(dat) {
  quantile(dat, probs = 0.50, na.rm = T)
}
