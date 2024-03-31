

#' Return level plots of generalized Pareto distribution (GPD)
#'
#' @param N  number of year for the return level one is seeking for
#' @param scale scale parameter of the GPD
#' @param shape scale parameter of the GPD
#' @param thr thresholds used to fit the GPD
#' @param m  number of days in a year
#' @param zeta pr(Y > thr)
#' @theta extremal index (1/theta corresponds to the cluster size)
#'
#' @return return level of the
#' @export
#'
#' @examples
return_l_gpd <- function(N, scale, shape, thr, m, zeta, theta) {
   if (abs(shape) > 1e-2) {
      z_N <- thr + (scale / shape) * (((N * m * zeta * theta) ^ shape) - 1)
   } else {
      z_N <- thr + scale * log(N * m * zeta * theta)
   }
   return(z_N)
}

#' Simulate and provide confidence interval of the return levels
#'
#' @param model.thr a list that define the models for thresholds to fit asymmetric Laplace distribution
#' @param model.gpd a list that define the models for GPD
#' @param zeta  1-zeta provide the quantiles ate which the thresholds are estimated using asymmetric Laplace distribution
#' @param data_train training data
#' @param data_test testing data
#' @param m  number of days in a year
#' @param zeta pr(Y > thr)
#' @theta extremal index (1/theta corresponds to the cluster size)
#' @N.year number of year for the return level one is seeking for
#' @nsim number of simulations for the estimation of the return level
#' @N.C2 return level for task C2
#' @return
#' @export
#'
#' @examples


fun_evgam_gpd <- function(fit_thr_ald,
                          model.gpd,
                          data_train,
                          data_test,
                          m,
                          zeta,
                          theta,
                          N.year,
                          nsim) {
   data_train$threshold <-
      predict(fit_thr_ald)$location ### the estimate of the thresholds based on the covariates
   data_train$excess <- data_train$Y - data_train$threshold
   is.na(data_train$excess[data_train$excess <= 0]) <-
      TRUE  ## put NAs to the non-threholds excesses which will be ignored when passing in evgam function
   ## fit of the GPD
   fit_gpd <- evgam::evgam(model.gpd, data_train, family = "gpd")
   
   #### simulation and prediction
   sim_ald <-
      simulate(fit_thr_ald,
               newdata = data_test,
               nsim = nsim,
               type = "response") ## simulation of the estimated thresholds using asymmetric Laplace
   sim_gpd <-
      simulate(fit_gpd,
               newdata = data_test,
               nsim = nsim,
               type = "response") ## simulation of the estimated scale and shape parameters of GPD
   
   return_levl <- matrix(NA, nrow = nsim, ncol = nrow(data_test))
   for (i in 1:nsim) {
      for (j in 1:nrow(data_test)) {
         return_levl[i, j] <- return_l_gpd(
            N = N.year,
            scale = sim_gpd[[1]][j, i],
            shape = sim_gpd[[2]][j, i],
            thr = sim_ald[[1]][j, i],
            m = m,
            zeta = zeta,
            theta = theta
         )
         
      }
   }
   CI_50 <-
      apply(
         return_levl,
         MARGIN = 2,
         FUN = quantile,
         probs = c(0.25, 0.75)
      )
   CI_95 <-
      apply(
         return_levl,
         MARGIN = 2,
         FUN = quantile,
         probs = c(0.025, 0.975)
      )
   Cond_quan_est <- apply(return_levl, MARGIN = 2, FUN = mean)
   C50_95_results <- cbind(
      Cond_quan_est = Cond_quan_est,
      lower_CI50 = CI_50[1, ],
      upper_CI50 = CI_50[2, ],
      lower_CI95 = CI_95[1, ],
      upper_CI95 = CI_95[2, ]
   )
   
   return(list("C50_95_results" = C50_95_results))
}
