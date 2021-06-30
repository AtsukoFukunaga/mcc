# Example of constructing regional MCC using pco scores after fitting covariates 

library(vegan)
library(tidyverse)
library(ggplot2)

source("res_pco_scores.R")
source("regional_mcc.R")

# load data (e.g. data file containing a data frame of fish count and a data frame for factors)
load("rda_files/sample_data_2.rda") # load "fish" and "factor" data frames

# get Bray-Curtis dissimilarity
fish_sqrt <- sqrt(fish_df)   # sqrt transform fish count, not always necessary to transform
fish_bray <- as.matrix(vegdist(fish_sqrt, "bray")) * 100  # a matrix of Bray-Curtis dissimilarity scaled from 0 to 100

# obtain pco scores
pco_data <- get_res_pco_scores(fish_bray, factor_df, factor_names = c("zone_depth")) 

# pco_df already has "region" column

# obtain dt and dbt for regional mcc with b = 2 (first 2 years as baseline data)
mcc_data <- regional_mcc(pco_data, b = 2)

# obtain bootstrap dt and dbt with 1000 bootstrap resampling
mcc_boot_data <- mcc_bootstrap(pco_data, mcc_data, b = 2, nboot = 1000)  # this will take some time to run

# dt data for plotting
dt <- data.frame(year = as.character(rep(c(2009:2016), 3)), region = rep(c("North", "Mid", "South"), each = 8))
dt_data <- mcc_boot_data$dt_data %>%
  rowwise() %>%
  mutate(dt95 = qt95(c_across(V1:V1000)), dt50 = qt50(c_across(V1:V1000))) %>%
  select(year, region, dt, dt95, dt50) %>% 
  right_join(dt)
dt95 <- dt_data %>%
  group_by(year) %>%
  summarise(dt95 = mean(dt95))
dt50 <- dt_data %>%
  group_by(year) %>%
  summarise(dt50 = mean(dt50))

# dbt data for plotting
dbt <- data.frame(year = as.character(rep(c(2009:2016), 3)), region = rep(c("North", "Mid", "South"), each = 8))
dbt_data <- mcc_boot_data$dbt_data %>%
  rowwise() %>%
  mutate(dbt95 = qt95(c_across(V1:V1000)), dbt50 = qt50(c_across(V1:V1000))) %>%
  select(year, region, dbt, dbt95, dbt50) %>%
  right_join(dbt)
dbt95 <- mean(dbt_data$dbt95, na.rm = TRUE)
dbt50 <- mean(dbt_data$dbt50, na.rm = TRUE)

ggplot(dt_data, aes(x = as.numeric(year), y = dt)) + 
  geom_point(aes(shape = region)) + 
  geom_line(aes(group = region)) + 
  xlab("Year") + 
  ylab(expression(italic(d[t]))) +
  ylim(0, 35) + 
  geom_line(data = dt95, aes(x = as.numeric(year), y = dt95), linetype = "twodash") + 
  geom_line(data = dt50, aes(x = as.numeric(year), y = dt50), linetype = "dashed")

ggplot(dbt_data, aes(x = as.numeric(year), y = dbt)) + 
  geom_point(aes(shape = region)) + 
  geom_line(aes(group = region)) + 
  xlab("Year") + 
  ylab(expression(italic({d^b}[t]))) + 
  ylim(0, 35) + 
  geom_hline(yintercept = dbt95, linetype = "twodash") + 
  geom_hline(yintercept = dbt50, linetype = "dashed")
  
