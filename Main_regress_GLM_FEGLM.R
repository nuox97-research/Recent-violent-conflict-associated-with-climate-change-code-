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

############################################################
# ==========================================================
# 1. GLM
# ==========================================================
############################################################
glm_origin <- feglm(
  binary_conflict ~ nicheloss_comb_lognorm,
  family = binomial(link = "logit"),
  data = Rawdata_model
)
summary(glm_origin)

# Moran's I index
resid_glm <- residuals(glm_origin, type = "pearson")
moran_glmm <- moran.test( resid_glm, lw, zero.policy = TRUE)
print(moran_glmm)

# consider SAC
glm_conley <- feglm(
  binary_conflict ~ nicheloss_comb_lognorm,
  family = binomial(link = "logit"),
  vcov = conley(500) ~ Lat + Lon,
  data = Rawdata_model
)
summary(glm_conley)

############################################################
# ==========================================================
# 2. FEGLM, considering niche loss only
# ==========================================================
############################################################
model_fe_conley_nicheloss <- feglm(
  binary_conflict ~ nicheloss_comb_lognorm | country,
  
  family = binomial(link = "logit"),
  vcov = conley(500)~Lat + Lon,
  data = Rawdata_model
)
summary(model_fe_conley_nicheloss)

############################################################
# ==========================================================
# 3. FEGLM
#    niche + niche loss
# ==========================================================
############################################################
model_fe_climate <- feglm(
  binary_conflict ~ niche_comb_norm + I(niche_comb_norm^2) + nicheloss_comb_lognorm | country,
  
  family = binomial(link = "logit"),
  vcov = conley(500)~Lat + Lon,
  data = Rawdata_model
)
summary(model_fe_climate)

############################################################
# ==========================================================
# 4. FEGLM
#    full model
# ==========================================================
############################################################
model_fe_conley_full <- feglm(
  binary_conflict ~ niche_comb_norm + I(niche_comb_norm^2) + nicheloss_comb_lognorm +
    GDPpercap_lognorm + langdiv_1dlognorm + gini_norm | country,
  
  family = binomial(link = "logit"),
  vcov = conley(500)~Lat + Lon,
  data = Rawdata_model
)
summary(model_fe_conley_full)


