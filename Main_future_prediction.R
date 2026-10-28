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

Rawdata_model <- fread("E:/DATA/Conflict/Code_2nd_version_additivenicheloss/Conflict/Conflict100km_n12707.csv")

############################################################
# ==========================================================
# Country FE + Conley Spatial SE.
# nicheloss only
# ==========================================================
############################################################
model_fe_conley <- feglm(
  binary_conflict ~ nicheloss_comb_lognorm | country,
  
  family = binomial(link = "logit"),
  vcov = conley(500)~Lat + Lon,
  data = Rawdata_model
)
summary(model_fe_conley)

# get the fixed effect coefficient
country_FE <- fixef(model_fe_conley) %>% unlist()
countryname <- sub("^country\\.", "", names(country_FE))
country_FE <- data.frame(
  country = countryname,
  FE = country_FE %>% as.numeric()
)
fwrite(country_FE,'country_FE_binary.csv')

############################################################
# ==========================================================
# Prediction
# Nicheloss only
# ==========================================================
############################################################
nicheloss_future = fread('Nicheloss_warm_100km.csv')
nicheloss_current <- Rawdata_model[, .(nicheloss_comb_lognorm = nicheloss_comb_lognorm,country,Lon,Lat)]
merged_data <- merge(nicheloss_current, nicheloss_future, by = c("Lon", "Lat"), all = FALSE)
fwrite(merged_data,'nicheloss_current_future.csv')

current <- merged_data[, .(nicheloss_comb_lognorm = nicheloss_comb_lognorm, country = country.x)] 
pred_data_warm15 <- merged_data[, .(nicheloss_comb_lognorm = nicheloss_warm15, country = country.x)]
pred_data_warm20 <- merged_data[, .(nicheloss_comb_lognorm = nicheloss_warm20, country = country.x)]
pred_data_warm24 <- merged_data[, .(nicheloss_comb_lognorm = nicheloss_warm24, country = country.x)]
pred_data_warm27 <- merged_data[, .(nicheloss_comb_lognorm = nicheloss_warm27, country = country.x)]
pred_data_warm36 <- merged_data[, .(nicheloss_comb_lognorm = nicheloss_warm36,country = country.x)]

pred_prob_current <- predict(model_fe_conley,newdata = current,type = "response")
pred_prob_warm15 <- predict(model_fe_conley,newdata = pred_data_warm15,type = "response")
pred_prob_warm20 <- predict(model_fe_conley,newdata = pred_data_warm20,type = "response")
pred_prob_warm24 <- predict(model_fe_conley,newdata = pred_data_warm24,type = "response")
pred_prob_warm27 <- predict(model_fe_conley,newdata = pred_data_warm27,type = "response")
pred_prob_warm36 <- predict(model_fe_conley,newdata = pred_data_warm36,type = "response")

Conflict_prob_future_change <- data.frame(
  Lon = merged_data$Lon,
  Lat = merged_data$Lat,
  probchange_warm15_abs = pred_prob_warm15 - pred_prob_current,
  probchange_warm20_abs = pred_prob_warm20 - pred_prob_current,
  probchange_warm24_abs = pred_prob_warm24 - pred_prob_current,
  probchange_warm27_abs = pred_prob_warm27 - pred_prob_current,
  probchange_warm36_abs = pred_prob_warm36 - pred_prob_current,
  country = merged_data$country.x
)
Conflict_prob_future_change <- na.omit(Conflict_prob_future_change)
fwrite(Conflict_prob_future_change,'Conflict_prob_future_nichelosschange_gridall_FEGLM.csv')


