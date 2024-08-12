library(evgam)

setwd(this.path::here())
load("../outputs/Amaurot_imputed.Rdata") ## imputed Amaurot data

#----------------------------------
# ALD model for threshold selection

model_thr_ald <- list(Y ~ Season + CP + s(V1, bs = "cs") + s(V2, bs = "cs") + 
                        s(V3, bs = "cs") + s(V4, bs = "cs") + 
                        s(WindDirection , bs = "cs") + s(WindSpeed , bs = "cs") + 
                        s(Atmosphere, bs = "cs"),
                      ~ Season + CP + s(V1, bs = "cs") + s(V2, bs = "cs") + 
                        s(V3, bs = "cs") + s(V4, bs = "cs") + 
                        s(WindDirection , bs = "cs") + s(WindSpeed , bs = "cs") + 
                        s(Atmosphere, bs = "cs"))

zeta.choices <- seq(0.95, 0.99, 0.01)

thr.levels.lists <- lapply(zeta.choices, function(zeta){
  fit_thr_ald <- evgam(formula = model_thr_ald, 
                       data = Amaurot.imp.final, 
                       family = "ald", ald.args = list(tau = zeta))
  fit_thr_ald})

save(thr.levels.lists, file = "../outputs/threshold_levels_lists.Rdata")

#----------------------------------
# list of all possible models for threshold exceedances

model_gpd.list <- list()

model_gpd.list[[1]] <- list(excess ~ 1, ~ 1)
model_gpd.list[[2]] <- list(excess ~ Season, ~ Season)
model_gpd.list[[3]] <- list(excess ~ Season + CP, ~ Season + CP)
model_gpd.list[[4]] <- list(excess ~ Season + CP + V1 + V2 + V3 + V4 +
                              WindDirection + WindSpeed + Atmosphere, ~ 1)
model_gpd.list[[5]] <- list(excess ~ Season + CP + V1 + V2 + V3 + V4 +
                              WindDirection + WindSpeed + Atmosphere, ~ Season)
model_gpd.list[[6]] <- list(excess ~ Season + CP + V1 + V2 + V3 + V4 +
                              WindDirection + WindSpeed + Atmosphere, ~ CP)
model_gpd.list[[7]] <- list(excess ~ Season + CP + V1 + V2 + V3 + V4 +
                              WindDirection + WindSpeed + Atmosphere, ~ Season + CP)

save(model_gpd.list, model_thr_ald, file = "../outputs/C1_models.Rdata")