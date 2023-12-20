library(parallel)
library(evgam)

setwd(this.path::here())
load("../outputs/Amaurot_imputed.Rdata")
load("../outputs/AmaurotTestSet_CP_10K.Rdata")
load("../outputs/models.Rdata")
load("../outputs/threshold_levels_lists.Rdata")

#-------------------------------------------------------------------------------

return_l_gpd <- function(N, scale, shape, thr, m, zeta, theta){
  if(abs(shape) > 1e-2){
    z_N<- thr + (scale/shape) * (((N * m * zeta * theta)^shape) -1)
  } else {
    z_N<- thr + scale * log(N * m * zeta * theta)  }
  return(z_N)
}

fun_evgam_gpd<- function(fit_thr_ald, 
                         model.gpd, 
                         data_train, 
                         data_test, 
                         m, 
                         zeta, 
                         theta, 
                         N.year, 
                         nsim){
  
  data_train$threshold<- predict(fit_thr_ald)$location 
  data_train$excess<- data_train$Y - data_train$threshold
  is.na(data_train$excess[data_train$excess<=0]) <- TRUE  
  ## fit of the GPD
  fit_gpd<- evgam::evgam(model.gpd, data_train, family = "gpd")
  
  #### simulation and prediction
  sim_ald <- simulate(fit_thr_ald, newdata = data_test, nsim = nsim, type = "response")
  sim_gpd <- simulate(fit_gpd, newdata = data_test, nsim = nsim, type = "response")
  
  return_levl<- matrix(NA, nrow = nsim, ncol = nrow(data_test))
  for (i in 1:nsim) {
    for (j in 1:nrow(data_test)) {
      return_levl[i,j]<- return_l_gpd(N=N.year,
                                      scale=sim_gpd[[1]][j,i], 
                                      shape=sim_gpd[[2]][j,i], 
                                      thr=sim_ald[[1]][j,i], 
                                      m=m, 
                                      zeta=zeta, 
                                      theta=theta)
      
    }
  }
  CI_50 <- apply(return_levl, MARGIN = 2, FUN = quantile, probs=c(0.25, 0.75))
  CI_95 <- apply(return_levl, MARGIN = 2, FUN = quantile, probs=c(0.025, 0.975))
  Cond_quan_est <- apply(return_levl, MARGIN = 2, FUN = mean)
  C50_95_results <- cbind(Cond_quan_est = Cond_quan_est, 
                          lower_CI50 = CI_50[1,],
                          upper_CI50 = CI_50[2,],
                          lower_CI95 = CI_95[1,],
                          upper_CI95 = CI_95[2,])
  
  return(list("C50_95_results" = C50_95_results))
}

zeta <- seq(0.95, 0.99, 0.01)
probs <- 0.9999 ## probability of interest
m <- 300 # number of days in a year
theta <- 1
nsim <- 1e4 # number of simulation
N.year <- 1 / (m * (1 - probs))  # number of years for the return level 

#---------------------------

results_thr_GP_models_95 <- parallel::mclapply(1:length(model_gpd.list), function(j){
  fun_evgam_gpd(fit_thr_ald=thr.levels.lists[[1]], 
                model.gpd=model_gpd.list[[j]], 
                data_train=Amaurot.imp.final, 
                data_test=AmaurotTestSet, 
                m=m, 
                zeta=1-zeta[1], 
                theta=theta, 
                N.year=N.year, 
                nsim=nsim)}, mc.cores = 10)

results_thr_GP_models_96 <- parallel::mclapply(1:length(model_gpd.list), function(j){
  fun_evgam_gpd(fit_thr_ald=thr.levels.lists[[2]], 
                model.gpd=model_gpd.list[[j]], 
                data_train=Amaurot.imp.final, 
                data_test=AmaurotTestSet, 
                m=m, 
                zeta=1-zeta[2], 
                theta=theta, 
                N.year=N.year, 
                nsim=nsim)}, mc.cores = 10)

results_thr_GP_models_97 <- parallel::mclapply(1:length(model_gpd.list), function(j){
  fun_evgam_gpd(fit_thr_ald=thr.levels.lists[[3]], 
                model.gpd=model_gpd.list[[j]], 
                data_train=Amaurot.imp.final, 
                data_test=AmaurotTestSet, 
                m=m, 
                zeta=1-zeta[3], 
                theta=theta, 
                N.year=N.year, 
                nsim=nsim)}, mc.cores = 10)

results_thr_GP_models_98 <- parallel::mclapply(1:length(model_gpd.list), function(j){
  fun_evgam_gpd(fit_thr_ald=thr.levels.lists[[4]], 
                model.gpd=model_gpd.list[[j]], 
                data_train=Amaurot.imp.final, 
                data_test=AmaurotTestSet, 
                m=m, 
                zeta=1-zeta[4], 
                theta=theta, 
                N.year=N.year, 
                nsim=nsim)}, mc.cores = 10)

results_thr_GP_models_99 <- parallel::mclapply(1:length(model_gpd.list), function(j){
  fun_evgam_gpd(fit_thr_ald=thr.levels.lists[[5]], 
                model.gpd=model_gpd.list[[j]], 
                data_train=Amaurot.imp.final, 
                data_test=AmaurotTestSet, 
                m=m, 
                zeta=1-zeta[5], 
                theta=theta, 
                N.year=N.year, 
                nsim=nsim)}, mc.cores = 10)

save(results_thr_GP_models_95, results_thr_GP_models_96, 
     results_thr_GP_models_97, results_thr_GP_models_98,
     results_thr_GP_models_99, file = "../outputs/C1_results_compare.Rdata")

#-------------------------------------------------------------------------------

COV50.table <- cbind(sapply(1:length(model_gpd.list), function(model.no){
  mean((results_thr_GP_models_95[[model.no]]$C50_95_results[ , 2] < AmaurotTestSet$Qtruth) * 
         (results_thr_GP_models_95[[model.no]]$C50_95_results[ , 3] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_96[[model.no]]$C50_95_results[ , 2] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_96[[model.no]]$C50_95_results[ , 3] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_97[[model.no]]$C50_95_results[ , 2] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_97[[model.no]]$C50_95_results[ , 3] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_98[[model.no]]$C50_95_results[ , 2] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_98[[model.no]]$C50_95_results[ , 3] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_99[[model.no]]$C50_95_results[ , 2] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_99[[model.no]]$C50_95_results[ , 3] > AmaurotTestSet$Qtruth))}))

#-----------------------------
# Coverage table (95%)

COV95.table <- cbind(sapply(1:length(model_gpd.list), function(model.no){
  mean((results_thr_GP_models_95[[model.no]]$C50_95_results[ , 4] < AmaurotTestSet$Qtruth) * 
         (results_thr_GP_models_95[[model.no]]$C50_95_results[ , 5] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_96[[model.no]]$C50_95_results[ , 4] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_96[[model.no]]$C50_95_results[ , 5] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_97[[model.no]]$C50_95_results[ , 4] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_97[[model.no]]$C50_95_results[ , 5] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_98[[model.no]]$C50_95_results[ , 4] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_98[[model.no]]$C50_95_results[ , 5] > AmaurotTestSet$Qtruth))}),
  sapply(1:length(model_gpd.list), function(model.no){
    mean((results_thr_GP_models_99[[model.no]]$C50_95_results[ , 4] < AmaurotTestSet$Qtruth) * 
           (results_thr_GP_models_99[[model.no]]$C50_95_results[ , 5] > AmaurotTestSet$Qtruth))}))

colnames(COV50.table) <- seq(0.95, 0.99, 0.01)
colnames(COV95.table) <- seq(0.95, 0.99, 0.01)

#-------------------------------------------------------------------------------

COV50_tab <- apply(COV50.table, 1:2, function(x){sprintf(x, fmt = "%.2f")})
COV95_tab <- apply(COV95.table, 1:2, function(x){sprintf(x, fmt = "%.2f")})

final_tab <- sapply(1:ncol(COV50_tab), function(j){
  sapply(1:nrow(COV50_tab), function(i){
    paste(COV50_tab[i, j], COV95_tab[i, j], sep = ", ")})})

rownames(final_tab) <- as.character(1:7)
colnames(final_tab) <- c("$q_{0.95}$","$q_{0.96}$","$q_{0.97}$","$q_{0.98}$","$q_{0.99}$")
cat(knitr::kable(final_tab, 
                 format = "latex", 
                 escape = FALSE, 
                 booktabs = TRUE),
    file = "../tables/Table4.tex")