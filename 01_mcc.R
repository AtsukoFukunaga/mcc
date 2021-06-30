library(vegan)
library(tidyverse)
library(ggplot2)


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


### MCC and bootstrap example

# load data (e.g. data file containing a data frame of fish count and a data frame for factors)
load("rda_files/data.rda") # load "fish" and "factor" data frames

# get Bray-Curtis dissimilarity
fish_sqrt <- sqrt(fish)   # sqrt transform fish count, not always necessary to transform
fish_bray <- as.matrix(vegdist(fish_sqrt, "bray")) * 100  # a matrix of Bray-Curtis dissimilarity scaled from 0 to 100

# obtain pco scores
pco_data <- get_pco_scores(fish_bray, factor)  # input should be a matrix of dissimilarity and a dataframe of factors with "year" and "site" columns
save(pco_data, file = "rda_files/fish_pco.rda")  # save pco_data to create regional mcc later

# obtain dt and dbt for mcc with b = 2 (first 2 years as baseline data)
mcc_data <- mcc(pco_data, b = 2)

# obtain bootstrap dt and dbt with 1000 bootstrap resampling
mcc_boot_data <- mcc_bootstrap(pco_data, mcc_data, b = 2, nboot = 1000)  # this will take some time to run

save(mcc_boot_data, file = "rda_files/fish_mcc.rda")  # give an appropriate name and save

rm(list = ls()) # empty the environment

load("rda_files/fish_mcc.rda")  # re-load the mcc_boot_data with 1000 bootstrap re-sampling

# dt data for plotting
dt_data <- mcc_boot_data$dt_data %>%
  rowwise() %>%
  mutate(dt95 = qt95(c_across(V1:V1000)), dt50 = qt50(c_across(V1:V1000))) %>%
  select(site, year, dt, dt95, dt50)

# dbt data for plotting
dbt_data <- mcc_boot_data$dbt_data %>%
  rowwise() %>% 
  mutate(dbt95 = qt95(c_across(V1:V1000)), dbt50 = qt50(c_across(V1:V1000))) %>%
  select(site, year, dbt, dbt95, dbt50) %>%
  group_by(site) %>%
  mutate(dbt95 = mean(dbt95, na.rm = TRUE), dbt50 = mean(dbt50, na.rm = TRUE))

# if plotting sites grouped by region, average dt95/dbt95 and dt50/dbt50 from the sites within each region for regional qt95 and qt50
