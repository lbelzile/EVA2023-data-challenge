
#' Loss function for task C2
#'
#' @param qhat estimated quantiles
#' @param q true quantiles
#'
#' @return
#' @export
#'
#' @examples
lossfun <- function(qhat, q){
   loss <- function(qhat, q){
      mean(ifelse(0.99*q > qhat,
                  0.99*(0.99*q-qhat),
                  ifelse(1.01*q < qhat,
                         0.1*(qhat-1.01*q),
                         0)))
   }
   retval <- numeric(length = length(qhat))
   for(i in seq_along(retval)){
      retval[i] <- loss(qhat[i], q = q)
   }
   retval
}


#' Density function of the asymmetric Laplace distribution
#'
#' @param x vector of quantile
#' @param location vector of location parameters
#' @param scale vector of scale parameters
#' @param tau asymmetric parameter giving the quantile level of \code{location}
#' @param log logical; if \code{TRUE}, return the log density
#' @return a vector of (log)-density evaluated at x
dald <- function(x, location, scale, tau, log = FALSE){
   stopifnot(length(tau) == 1L,
             isTRUE(all(tau < 1, tau > 0, scale > 0)),
             length(location) %in% c(1L, length(x)),
             length(scale) %in% c(1L, length(x)))
   yc <- (x-location)/scale
   logdens <- log(tau) + log(1-tau) - log(scale) -yc*(tau-I(yc < 0))
   if(log){
      return(logdens)
   } else{
      return(exp(logdens))
   }
}
#' Distribution function of the asymmetric Laplace distribution
#'
#' @param q vector of quantiles
#' @inheritParams dald
#' @param lower.tail logical; if \code{TRUE}, returns the survival probability
#' @param log.p logical; if \code{TRUE}, returns the log probability
#' @return vector of distribution function evaluated at \code{q}
pald <- function(q, location, scale, tau, lower.tail = TRUE, log.p = FALSE){
   yc <- (q-location)/scale
   prob <- ifelse(yc < 0, tau*exp((1-tau)*yc), 1-(1-tau)*exp(-tau*yc))
  if(!lower.tail){
     prob <- 1-prob
  }
   if(log.p){
      return(log(prob))
   } else{
      return(prob)
   }
}


#' Bayesian bootstrap
#'
#' Given a data frame of inputs, create a copy
#' with nsize rows and sample with replacement
#' using Dirichlet(1) weights for each element
#' of groups.
#' @param data a data frame of observations
#' @param nsize number of rows of the output
#' @param group a list with covariate names.
#'
#' @note  All of the column names of \code{data} must be given in \code{groups}
#' @return A \code{nsize} by \code{ncol(data)} data frame
bayes_boot <- function(data, nsize,
                       groups = list(c("V1","V2"), c("Season","V3"), "V4",
                                     c("WindDirection","WindSpeed","CP"),
                                     "Atmosphere")){
 n <- nrow(data)
 newdata <- dplyr::slice_sample(data, n = nsize)
 stopifnot(isTRUE(all(colnames(data) %in% unlist(groups))))
 for(j in seq_along(groups)){
    col_ids <- which(groups[[j]] %in% colnames(newdata))
    prob <- revdbayes::rDir(n = 1, alpha = rep(1, n))[1,]
    obs <- sample.int(n = n, size = nsize, replace = TRUE, prob = prob)
    newdata[,col_ids] <- data[obs, col_ids]
}
   return(newdata)
}

#' Calculate loss function for nonstationary generalized Pareto models
#'
#' @param fit_thr_ald fitted EVGAM object for threshold exceedances
#' @param model.gpd a list that define the generalized Pareto distribution model
#' @param zeta  1-zeta provide the quantiles at which the thresholds are estimated using asymmetric Laplace distribution; i.e., pr(Y > thr)
#' @param data_train training data: a matrix of response and covariates; the first column is the response
#' @param newdata matrix of new covariate values
#' @param bayesboot logical; if \code{TRUE}, sample covariates using a Bayesian bootstrap
#' @param qvals values at which to estimate the loss
#' @param m number of days in a year.
#' @param nsim number of posterior draws for the estimation of the return level
#' @param ncov number of covariates to subsample
#' @param qvals quantile values at which to evaluate the loss function
#' @param Nyear number of year for the return level
#' @return
#' @export
fun_evgam_gpd <- function(fit_thr_ald,
                          model.gpd,
                          data_train,
                          newdata,
                          method = c("npboot","bayesboot","sweep"),
                          m = 300L,
                          zeta = 1 - fit_thr_ald$tau,
                          Nyear = 200L,
                          nsim = 1e4L,
                          ncov = 1e4L,
                          seed = 2023,
                          retlev.method = c("Eastoe.Tawn", "Youngman"),
                          fix_zeta = FALSE) {
   set.seed(seed)
   retlev.method <- match.arg(retlev.method)
   ### the estimate of the thresholds based on the covariates
   data_train$threshold <- predict(fit_thr_ald, newdata = data_train)$location
   data_train$excess <- data_train$Y - data_train$threshold
   ## put NAs to the non-threshold excesses
   data_train$excess <- with(
      data_train,
      ifelse(excess <= 0, NA, excess))
   ## fit GP distribution to threshold exceedances
   fit_gpd <- evgam::evgam(formula = model.gpd,
                           data = data_train,
                           family = "gpd")
   # Simulate nsim values from the posterior and compute associated
   # thresholds, scale and shape for each
   retlev <- as.numeric(replicate(n = ceiling(nsim/100), expr = {
      if(method == "bayesboot"){
         ndat <- bayes_boot(newdata, nsize = ncov)
      } else if(method == "npboot"){
         ndat <- dplyr::slice_sample(newdata, n = ncov)
      } else{
         # Sweep over all existing data
         ndat <- newdata
      }

   # Generate thresholds from fitted ALD model and
   # obtain approximate posterior pred of zeta_u

   nsimu <- 100
   if(fix_zeta){
      # Ignore probability of exceedance and keep the value constant
      pthresh <- 1 - zeta
      # Instead, use the location parameter to vary the threshold
      u <- simulate(
         fit_thr_ald,
         newdata = ndat,
         nsim = nsimu,
         type = "response")$location
   } else{
   # Sample from approximate posterior of the ALD parameters
   # via a Gaussian approximation on the transformed scale
   u <- predict(fit_thr_ald, newdata = ndat)[,'location']
   post_ald <- simulate(
         fit_thr_ald,
         newdata = ndat,
         nsim = nsimu,
         type = "response")
   # Compute the probability of exceeding u
   # treating the latter as a fixed quantity
   pthresh <- sapply(seq_len(nsimu), function(i){
      pald(q = u,
        location = post_ald$location[,i],
        scale = post_ald$scale[,i],
        tau = 1 - zeta)})
   }
   # Simulate generalized Pareto parameters and compute
   # resulting scale and shape for each combination of covariates
   gp_sim <- simulate(
       fit_gpd,
       newdata = ndat,
       nsim = nsimu,
       type = "response")
   # Given matrix of inputs, compute the return levels associated
   if(retlev.method == "Youngman"){
   retlev <- evgam::qev(
         p = 1 - (1 / Nyear),
         loc = u,
         scale = gp_sim[[1]],
         shape = gp_sim[[2]],
         m = m,
         family = "gpd",
         tau = pthresh
   )
   } else if(retlev.method == "Eastoe.Tawn"){
    retlev <- sapply(
       seq_len(nsimu),
       function(j){
         suppressWarnings(uniroot(f = function(x){
            log(mean((1-pthresh[,j])*
               revdbayes::pgp(
                    q = x,
                    loc = u,
                    scale = gp_sim[[1]][,j],
                    shape = gp_sim[[2]][,j],
                    log.p = FALSE,
                    lower.tail = FALSE))) +
               log(Nyear*m)
                 }, interval = c(50, 500))$root)})
   }
     }, simplify = TRUE))
   loss <- optimize(f = function(q){
      lossfun(qhat = q, q = retlev)
      }, interval = c(0.9, 1.4)*range(retlev), tol = 1e-3)$minimum

   ## Plot curves
   # xgrid <- seq(0.9*min(retlev), 1.1*max(retlev), length.out = 1001L)
   # lossfn <- sapply(xgrid,
   #    function(q){
   #       lossfun(qhat = q, q = retlev)
   #    })
   # plot(xgrid, lossfn, type = "l")
   # return(list(quantile = xgrid, lossfn, minloss = loss))
   return(loss)
}

#' Function to compute return levels
#'
#' This function was used as a sanity check to compare the results of qev
#' and the model output for the generalized Pareto and the T-year maximum
#' formulation. The latter works only if the thresholds are not too high.
#'
#' This function returns the same quantity as \code{qev} in package \code{evgam}.
#' @param loc matrix of threshold parameters
#' @param scale matrix of generalized Pareto scale
#' @param shape matrix of generalized Pareto shape
#' @param zeta matrix or scalar probability of exceedance above the threshold
#' @param Nyear number of years for the return levels
#' @param m number of observations per year
#' @param method string; whether to use the \code{GPD} to compute return levels or the \code{Tmax} distribution approximation
#' @param interval vector of length 2 containing the bounds for the root search.
#' @return a vector of \code{ncol(loc)} return levels.
retlev_fn <- function(
      loc,
      scale,
      shape,
      zeta,
      Nyear,
      m,
      method = c("GPD", "Tmax"),
      interval = c(150, 300)){
   method <- match.arg(method)
   stopifnot(is.matrix(loc), is.matrix(scale), is.matrix(shape),
             nrow(scale) == nrow(shape), nrow(loc) == nrow(shape),
             ncol(scale) == ncol(shape), ncol(loc) == ncol(shape))
   nc <- ncol(scale)
   retlev <- numeric(length = nc)
   for(j in 1:ncol(loc)){
      if(method == "Tmax"){
         retlev[j] <- uniroot(f = function(x){
            Nyear*zeta*m*log(mean(revdbayes::pgp(
               q = x, loc = loc[,j], scale = scale[,j],
               shape = shape[,j]))) - log(0.368)},
            interval = interval)$root
      } else if(method == "GPD"){
         retlev[j] <- uniroot(f = function(x){
            m*log(mean(1-zeta*revdbayes::pgp(
               q = x, loc = loc[,j], scale = scale[,j],
               shape = shape[,j], lower.tail = FALSE))) - log(1-1/Nyear)},
            interval = interval)$root
      }
   }
   return(retlev)
}
