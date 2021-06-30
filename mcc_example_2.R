# Example of constructing regional MCC using get_pco_scores() 

library(vegan)
library(ggplot2)

source("pco_scores.R")
source("regional_mcc.R")

# load data (e.g. data file containing a data frame of fish count and a data frame for factors)
load("rda_files/sample_data_1.rda") # load "fish" and "factor" data frames

# get Bray-Curtis dissimilarity
fish_sqrt <- sqrt(fish_df)   # sqrt transform fish count, not always necessary to transform
fish_bray <- as.matrix(vegdist(fish_sqrt, "bray")) * 100  # a matrix of Bray-Curtis dissimilarity scaled from 0 to 100

# obtain pco scores
pco_data <- get_pco_scores(fish_bray, factor_df)  # input should be a matrix of dissimilarity and a dataframe of factors with "year" and "site" columns

# add "region" column to pco_df if it does not exist
# e.g.
# "F01" "F02" "F03" "F04" "F05" in "north" region
# "F06" "F07" "F08" "F09" "F10" in "middle" region
# "F11" "F12" "F13" "F14" "F15" in "south" region
pco_data$pco_df$region <- ""
pco_data$pco_df$region[pco_data$pco_df$site %in% c("F01", "F02", "F03", "F04", "F05")] <- "north"
pco_data$pco_df$region[pco_data$pco_df$site %in% c("F06", "F07", "F08", "F09", "F10")] <- "middle"
pco_data$pco_df$region[pco_data$pco_df$site %in% c("F11", "F12", "F13", "F14", "F15")] <- "south"

# obtain dt and dbt for regional mcc with b = 2 (first 2 years as baseline data)
mcc_data <- regional_mcc(pco_data, b = 2)

# obtain bootstrap dt and dbt with 1000 bootstrap resampling
mcc_boot_data <- mcc_bootstrap(pco_data, mcc_data, b = 2, nboot = 1000)  # this will take some time to run

# dt data for plotting
dt_data <- mcc_boot_data$dt_data %>%
  rowwise() %>%
  mutate(dt95 = qt95(c_across(V1:V1000)), dt50 = qt50(c_across(V1:V1000))) %>%
  select(year, region, dt, dt95, dt50)

# dbt data for plotting
dbt_data <- mcc_boot_data$dbt_data %>%
  rowwise() %>%
  mutate(dbt95 = qt95(c_across(V1:V1000)), dbt50 = qt50(c_across(V1:V1000))) %>%
  select(year, region, dbt, dbt95, dbt50) %>%
  group_by(region) %>%
  mutate(dbt95 = mean(dbt95, na.rm = TRUE), dbt50 = mean(dbt50, na.rm = TRUE))

ggplot(dt_data, aes(x = as.numeric(year), y = dt)) + 
  geom_point() + 
  geom_line() + 
  xlab("Year") + 
  ylab(expression(italic(d[t]))) + 
  geom_line(data = dt_data, aes(x = as.numeric(year), y = dt95), linetype = "twodash") + 
  geom_line(data = dt_data, aes(x = as.numeric(year), y = dt50), linetype = "dashed") + 
  facet_wrap(~region, ncol = 1)

ggplot(dbt_data, aes(x = as.numeric(year), y = dbt)) + 
  geom_point() + 
  geom_line() + 
  xlab("Year") + 
  ylab(expression(italic({d^b}[t]))) + 
  geom_line(data = dbt_data, aes(x = as.numeric(year), y = dbt95), linetype = "twodash") + 
  geom_line(data = dbt_data, aes(x = as.numeric(year), y = dbt50), linetype = "dashed") + 
  facet_wrap(~region, ncol = 1)
