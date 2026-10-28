rm(list = ls())

# Here, we illustrate the approach using all binary conflicts
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

# group by country
Rawdata_model <- Rawdata_model %>%
  group_by(country) %>%
  mutate(
    nicheloss_between = mean(nicheloss_comb_lognorm, na.rm = TRUE),
    nicheloss_within = nicheloss_comb_lognorm - nicheloss_between
  ) %>%
  ungroup() %>%
  data.frame()

############################################################
# between-within decomposition 
############################################################
glm_conley <- feglm(
  binary_conflict ~ nicheloss_within + nicheloss_between,
  family = binomial(link = "logit"),
  vcov = conley(500) ~ Lat + Lon,
  data = Rawdata_model
)
summary(glm_conley)


