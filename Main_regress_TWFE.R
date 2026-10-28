rm(list = ls())

############################################################
# Packages
############################################################

library(data.table)
library(dplyr)
library(sf)

library(glmmTMB)
library(fixest)

library(spdep)
library(parallel)

set.seed(123)

############################################################
# Data
############################################################

setwd(getwd())

Allgrid <- fread("Conflict100km_n12707_nichelossorigin.csv")
Allgrid <- na.omit(Allgrid)
Allgrid$Pop[Allgrid$Pop == 0] <- 0.01

Conflict_event <- fread("Conflict_event_grid_nichelossorigin.csv")
Conflict_event$Pop[Conflict_event$Pop == 0] <- 0.01

###################################################################

################# discrete-time binary panel model ################

###################################################################
yearstart = 1990
yearall <- 1990:2020
grid  <- unique(Allgrid$gridid)

panel <- expand.grid(gridid = grid, year = yearall)

panel <- merge(panel, Allgrid, by = "gridid", all.x = TRUE)
Conflict_count <- Conflict_event[
  ,
  .(n_event = .N),
  by = .(gridid = id_2, year)
]
panel <- merge(panel,Conflict_count,by = c("gridid", "year"), all.x = TRUE)
panel$n_event[is.na(panel$n_event)] <- 0

# linearization
panel$bievent <- panel$n_event>0 %>% as.numeric()
panel$nicheloss_t <- scale(log((panel$nicheloss_comb_origin * (panel$year - yearstart + 1) / length(yearall))+1)) %>% as.vector()
panel$conpercaparea <- scale(log(panel$n_event/panel$Pop/panel$pixelArea*10^10+1)) %>% as.vector()
panel$condensity <- scale(log(panel$n_event/panel$pixelArea*10^5+1)) %>% as.vector()

# for binary conflict
model_binary <- feglm(
  bievent ~ nicheloss_t | country + year,
  family = binomial("logit"),
  vcov = conley(500)~Lat + Lon,
  data = panel
)
summary(model_binary)

# for conflict density
model_density <- feglm(
  condensity ~ nicheloss_t| country + year,
  family = gaussian(),
  vcov = conley(500)~Lat + Lon,
  data = panel
)
summary(model_density)

# for conflict density per capita
model <- feglm(
  conpercaparea ~ nicheloss_t | country + year,
  family = gaussian(),
  vcov = conley(500)~Lat + Lon,
  data = panel
)
summary(model)
