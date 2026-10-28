rm(list = ls())

# Here, we illustrate the approach using binary conflicts
# for conflict density (or conflict density per capita), replace the response variable with conflictden_lognorm (or conflict_percaparea_lognorm)
# and switch the model family from binomial to gaussian

############################################################
# Packages
############################################################

library(data.table)
library(dplyr)
library(sf)
library(glmmTMB)
library(fixest)
library(spdep)

############################################################
# Data
############################################################

setwd(getwd())

Rawdata_model <- fread("Conflict100km_n12707.csv")

############################################################
# model
############################################################
model_comp <- feglm(
  binary_conflict ~ nicheloss_comb_lognorm + GDPpercap_lognorm | country,
  data = Rawdata_model,
  vcov = conley(500) ~ Lat + Lon,
  family = binomial(link="logit")
)
summary(model_comp)

beta_niche <- coef(model_comp)["nicheloss_comb_lognorm"]
beta_GDP   <- coef(model_comp)["GDPpercap_lognorm"]

# GDP SD increase required to offset 1 SD increase in niche loss
comp_ratio <- -beta_niche / beta_GDP
comp_ratio

############################################################
# Delta method confidence interval
############################################################
V <- vcov(model_comp)
var_ratio <- 
  (1/beta_GDP)^2 * V["nicheloss_comb_lognorm","nicheloss_comb_lognorm"] +
  (beta_niche/beta_GDP^2)^2 * V["GDPpercap_lognorm","GDPpercap_lognorm"] -
  2*(1/beta_GDP)*(beta_niche/beta_GDP^2)*
  V["nicheloss_comb_lognorm","GDPpercap_lognorm"]

se_ratio <- sqrt(var_ratio)

ratio_lower <- comp_ratio - 1.96*se_ratio
ratio_upper <- comp_ratio + 1.96*se_ratio

data.frame(
  compensation_SD = comp_ratio,
  lower95 = ratio_lower,
  upper95 = ratio_upper
)

############################################################
# Convert SD compensation into GDP percentage change
############################################################
sd_logGDP <- 1.0296
delta_logGDP <- comp_ratio * sd_logGDP
GDP_percent_change <- 
  (exp(delta_logGDP)-1)*100
GDP_percent_change

############################################################
# Scenario compensation
############################################################
nicheloss_current_future = fread('nicheloss_current_future.csv')

warming <- data.frame(
  deltanicheloss15=
    nicheloss_current_future$nicheloss_warm15 -
    nicheloss_current_future$nicheloss_comb_lognorm,
  deltanicheloss20=
    nicheloss_current_future$nicheloss_warm20 -
    nicheloss_current_future$nicheloss_comb_lognorm,
  deltanicheloss27=
    nicheloss_current_future$nicheloss_warm27 -
    nicheloss_current_future$nicheloss_comb_lognorm,
  deltanicheloss36=
    nicheloss_current_future$nicheloss_warm36 -
    nicheloss_current_future$nicheloss_comb_lognorm
)

GDP_SD_required <- data.frame(
  GDP_SD_required_15 = -warming$deltanicheloss15 * beta_niche / beta_GDP,
  GDP_SD_required_20 = -warming$deltanicheloss20 * beta_niche / beta_GDP,
  GDP_SD_required_27 = -warming$deltanicheloss27 * beta_niche / beta_GDP,
  GDP_SD_required_36 = -warming$deltanicheloss36 * beta_niche / beta_GDP
)

GDP_percent_required<- data.frame(
  Lon = nicheloss_current_future$Lon,
  Lat = nicheloss_current_future$Lat,
  GDP_percent_required_15 = (exp(GDP_SD_required$GDP_SD_required_15 * sd_logGDP)-1)*100,
  GDP_percent_required_20 = (exp(GDP_SD_required$GDP_SD_required_20 * sd_logGDP)-1)*100,
  GDP_percent_required_27 = (exp(GDP_SD_required$GDP_SD_required_27 * sd_logGDP)-1)*100,
  GDP_percent_required_36 = (exp(GDP_SD_required$GDP_SD_required_36 * sd_logGDP)-1)*100
)
fwrite(GDP_percent_required,'GDP_percent_required_binary.csv')

