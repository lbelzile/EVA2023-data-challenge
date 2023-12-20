setwd(this.path::here())
library(evgam)
load("../outputs/threshold_levels_lists.Rdata")
load("../outputs/models_and_data_imputed.Rdata")
load("../outputs/Amaurot_imputed.Rdata")
Amaurot <- read.csv("../data/Amaurot.csv")

# These files contain the imputation and fitted models
# for the single imputation - corresponding to our submission
load("../outputs/train_test_data_imputed.Rdata")
load("../outputs/fitted-thr-models_and_GP-models.Rdata")
source("Task2-functions.R")


mods_id <- 1:7 #c(1,18, 42:46)
thresh_id <- 1:4
nrep <- 10L
# Try evaluating only on complete data
# We need to load the imputed data because it includes the 'Changepoint'
# variable for wind direction regime.
newdata <- Amaurot.train[
   apply(Amaurot, 1, function(x){!isTRUE(any(is.na(x)))}),]
if(!"outputs/C2_results_finalsubmission.Rdata" %in% list.files("..", recursive = TRUE)){
task2_results <-  array(dim = c(length(thresh_id), length(mods_id), nrep))
for(i in seq_along(thresh_id)){
   for(j in seq_along(mods_id)){
      for(k in seq_len(nrep)){
         ### loop over the GP models
         loss <-
            try(fun_evgam_gpd(
               fit_thr_ald = list_fit_thr_ald[[i]],
               ### output of the fitted thresholds models i
               model.gpd = model_gpd.list[[mods_id[j]]],
               ### the formula for the GPD models j
               data_train = Amaurot.train,
               newdata = newdata,
               method = "npboot",
               fix_zeta = FALSE,
               seed = 2023 + (k - 1L)
            ))
         if(!inherits(loss, "try-error")){
            task2_results[i,j,k] <- loss
         }
         print(paste0("Threshold ", i,
                      ", model ", j,
                      ", repetition", k,
                      ": ", Sys.time()))
      }
      save(task2_results, file = "../outputs/C2_results_finalsubmission.Rdata")
   }
}
} else{
   load("../outputs/C2__finalsubmission.Rdata")
}
retlev <- apply(task2_results, 1:2, mean)
retlev_sd <- apply(task2_results, 1:2, sd)/sqrt(9)
tab <- matrix("", nrow = ncol(retlev), ncol = nrow(retlev))
for(i in seq_len(nrow(task2_results))){
   for(j in seq_len(ncol(task2_results))){
     tab[j,i] <-  paste0(
        sprintf(round(retlev[i,j],1), fmt = "%.1f"), " (",
             sprintf(retlev_sd[i,j], fmt = "%.2f"),")")
   }
}
rownames(tab) <- 1:nrow(tab)
cat(
   knitr::kable(tab,
                booktabs = TRUE,
                format = "latex",
                col.names = paste0("$q_{", c(0.95, 0.96, 0.97, 0.98), "}$"),
                linesep = "",
                row.names = TRUE,
                escape = FALSE,
                align = c('rrrr')),
   file = "../tables/Table5.tex")

# In the final submission, we picked a single model
# We fixed zeta_u to the MLE
# We computed the marginal by sweeping over all existing observations
# (including imputed values) rather than resampling with replacement from the complete cases
if(!"outputs/C2_results_submission.Rdata" %in% list.files("..", recursive = TRUE)){
   results_thr_GP_models_sub <- matrix(nrow = length(thresh_id), ncol = length(mods_id))
   for(i in seq_along(thresh_id)){
      for(j in seq_along(mods_id)){
         ### loop over the GP models
         loss <-
            try(fun_evgam_gpd(
               fit_thr_ald = thr.levels.lists[[i]],
               ### output of the fitted thresholds models i
               model.gpd = model_gpd.list[[mods_id[j]]],
               ### the formula for the GPD models j
               data_train = Amaurot.train,
               newdata = Amaurot.train,
               method = "sweep",
               fix_zeta = TRUE
            ))
         if(!inherits(loss, "try-error")){
            results_thr_GP_models_sub[i,j] <- loss
         }
         print(paste0("Threshold ", i, ", model ", j, ": ", Sys.time()))
      }
   }
   save(results_thr_GP_models_sub, file = "../outputs/C2_results_submission.Rdata")
} else{
   load("../outputs/C2_results_submission.Rdata")
}
retlev <- matrix(unlist(unlist(results_thr_GP_models_sub)), nrow = length(mods_id))
tab <- round(retlev, 1)
rownames(tab) <- 1:nrow(tab)
cat(
   knitr::kable(tab,
                booktabs = TRUE,
                format = "latex",
                col.names = paste0("$q_{", c(0.95, 0.96, 0.97, 0.98), "}$"),
                linesep = "",
                row.names = TRUE,
                escape = FALSE),
   file = "../tables/Task2-submission.tex")
