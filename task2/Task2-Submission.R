setwd(this.path::here())

load("../outputs/C2_results_compare.Rdata")

tau <- c(0.95, 0.955, 0.96, 0.965, 0.97, 0.975, 0.98)
qu <- seq(150, 250, by = 0.01)


loss.thr.models <-
   loglik <- aic <- bic <- matrix(NA, nrow = 39, ncol = length(tau))
for (i in 1:7) {
   ### loop over possible thresholds
   for (j in 1:39) {
      ### loop over generalized Pareto models
      retlev_gp <-
         results_thr_GP_models[[i]][["results_GP_models"]][[j]]$r_200_gpd_sample
      loss <- results_thr_GP_models[[i]][["lossC2"]][[j]]
      loss.thr.models[j, i] <- qu[which.min(loss)]
      loglik[j, i] <-
         results_thr_GP_models[[i]][["results_GP_models"]][[j]]$log_lik
      aic[j, i] <-
         results_thr_GP_models[[i]][["results_GP_models"]][[j]]$aic_vale
      bic[j, i] <-
         results_thr_GP_models[[i]][["results_GP_models"]][[j]]$bic_val
   }
}

### model comparisons
rownames(loss.thr.models) <- paste0("Model-", 1:39)
colnames(loss.thr.models) <- paste0("tau", tau)
write.csv(loss.thr.models, file = "C2_usingLossFunctions.csv")

### results to submit
AnswerC2 <- qu[which.min(results_thr_GP_models[[7]][[2]][[39]])]
AnswerC2
