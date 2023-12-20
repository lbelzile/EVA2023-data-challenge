#' Empirical estimation of the tail probability
#'
#' @param data \code{n} by 3 data matrix for the Heffernan-Tawn model, in Laplace margins
#' @param levels above which to compute the levels. Default to values on Gumbel margins
#' @param prob string; the probability to estimate, either \code{p1} or \code{p2}.
#' @return estimated probability
#' @export
emp_C3 <- function(data, levs = c(6, 7, -log(log(2))), prob = c("p1","p2")) {
   prob <- match.arg(prob)
   Y1 <- data[, 1]
   Y2 <- data[, 2]
   Y3 <- data[, 3]
   if(prob == "p1"){
      y <- levs[1]
      mean((Y1 > y) & (Y2 > y) & (Y3 > y))
   } else{
      v <- levs[2]
      m <- levs[3]
      mean((Y1 > v) & (Y2 > v) & (Y3 < m))
   }
}

#' Conditional multivariate extreme value model of Heffernan and Tawn (2004) for Task 3
#'
#' @param data data matrix of dimension \code{n} by 3
#' @param mqu marginal quantile probability
#' @param cond.var index of the conditioning variable
#' @param dqu dependence quantile. Used to specify the marginal quantile at which to threshold the conditioning variable data when estimating the dependence parameters
#' @param nsim number of simulated example from the fitted model
#'
#' @return simulated observation from the fitted model and estimated required probabilities
#' @export
#'
#' @examples
probs_texmex_HT <-
   function(data, threshold, nsim = 1e4L, mqu = 0.9) {
      yp <- exp(-exp(-6))
      y <- qlaplace(yp)
      vp <- exp(-exp(-7))
      v <- qlaplace(vp)
      mp <- 0.5
      m <- 0
      dqu <- threshold
      mqu <- rep(mqu, length.out = ncol(data))
      ## transform data to Laplace margins
      fit1 <- texmex::mex(
            data,
            mqu = mqu,
            dqu = dqu,
            which = 1,
            penalty = "none"
         )
      fit2 <- texmex::mex(
         data,
         mqu = mqu,
         dqu = dqu,
         which = 2,
         penalty = "none"
      )
      fit3 <- texmex::mex(
         data,
         mqu = mqu,
         dqu = dqu,
         which = 3,
         penalty = "none"
      )
      plev <- 0.97
      ## simulation
      pred1_1 <- predict(fit1,
                 which = 1,
                 pqu = plev,
                 nsim = nsim)
      pred1_2 <- predict(fit2,
                 which = 2,
                 pqu = plev,
                 nsim = nsim)
      pred1_3 <- predict(fit3,
                         which = 3,
                         pqu = plev,
                         nsim = nsim)
      ## simulation based on point estimates
      ## summary of predicted values
      #sum.predict<- summary(pred, probs=c(0.01,0.2,0.5,0.8,0.99))
      probs.est_1 <- (1-plev) * mean(c(mean(rowSums(pred1_1$data$simulated > y) == 3),
                                       mean(rowSums(pred1_2$data$simulated > y) == 3),
                                       mean(rowSums(pred1_3$data$simulated > y) == 3)))
      pred2_1 <- predict(fit1,
                         which = 1,
                         pqu = plev,
                         nsim = nsim)
      pred2_2 <- predict(fit2,
                         which = 2,
                         pqu = plev,
                         nsim = nsim)
      prob2 <- function(data){
         mean(apply(data, 1, function(x){isTRUE(all(x[1] > v & x[2] > v & x[3] < m))}))
      }
      probs.est_2 <- (1-plev)* mean(c(
         prob2(pred2_1$data$simulated),
         prob2(pred2_2$data$simulated)))
      c(probs.est_1, probs.est_2)
}


#' Fit the conditional extremes model
#'
#' @param data data on standard Laplace margins
#' @param cond.var integer giving the index of the conditioning variable
#' @param thresh threshold on the probability scale
#' @return a list with elements
#' \itemize{
#' \item \code{cond.var} index of the conditioning variable
#' \item \code{threshold} threshold on the standard Laplace scale
#' \item \code{alpha} vector of parameters
#' \item \code{beta} vector of parameters
#' \item \code{resid} matrix of residuals from the fit
#' }
fit_HT <- function(data, cond.var, thresh){
   u <- qlaplace(thresh)
   stopifnot(length(cond.var) == 1L, cond.var <= ncol(data), cond.var >= 1L)
   X <- data[data[, cond.var] > u, -cond.var, drop = FALSE]
   X0 <- data[data[, cond.var] > u, cond.var]
   p <- ncol(data) - 1L
   ### objective function to be minimized
   fht_pll <- function(par, X, X0) {
      # Define parameters
      p <- ncol(as.matrix(X))
      alpha <- par[1:p]
      beta <- par[1:p + p]
      loc <- par[1:p + 2*p]
      scale <- exp(par[1:p + 3*p])
      obj <- 0
      for (j in seq_len(ncol(X))) {
         obj <- obj + sum(dnorm(
            X[, j],
            mean = alpha[j] * X0 + scale[j] * X0^beta[j] * loc[j],
            sd = scale[j] * X0^beta[j],
            log = TRUE
         ))
      }
      if(!is.finite(obj)){
         obj <- -1e10
      }
      return(-obj)
   }
   p <- ncol(as.matrix(X))
   #### optimization to find the parameter estimates
   opt <- Rsolnp::solnp(
      pars = c(rep(0.2, p), rep(0.1, p), rep(0, 2*p)),
      fun = fht_pll,
      LB = c(rep(-1, p), rep(-1e3, 3*p)),
      UB = c(rep(1, 2*p), rep(1e3, 2*p)),
      control = list(trace = FALSE),
      X = X, X0 = X0
   )
   if(opt$convergence != 0){
      warning("Optimization did not converge.")
   }
   alpha <- opt$pars[1:p]
   beta <- opt$pars[1:p + p]
   loc <- opt$pars[1:p + 2*p]
   scale <- exp(opt$pars[1:p + 3*p])
   Z <- matrix(0, nrow = nrow(X), ncol = ncol(X))
   for(j in seq_len(p)){
      Z[,j] <- (X[,j] - (alpha[j] * X0 + scale[j] * X0^beta[j] * loc[j]))/(scale[j] * X0^beta[j])
      Z[,j] <- loc[j] + Z[,j]*scale[j]
   }
   return(
      list(
         cond.var = cond.var,
         u = u,
         alpha = alpha,
         beta = beta,
         resid =  Z)
   )
}

#' Quantile function of the Laplace distribution
#' 
#' @param p a vector or matrix of probabilities
#' @return a vector or matrix of quantiles
qlaplace <- function(p) {
   -sign(p-0.5) * log(2 * pmin(p, 1-p))
}

#' Distribution function of the Gumbel distribution
#' 
#' @param q a vector or matrix of quantiles
#' @return a vector or matrix of probabilities
pgumbel <- function(q){
   exp(-exp(-q))
}

#' Map Gumbel data to standard Laplace scale
#' 
#' @param data a vector or matrix of quantiles in standard Gumbel margins
#' @return a vector or matrix of quantiles  in standard Laplace margins
gumbel2laplace <- function(data){
   qlaplace(pgumbel(data))
}

#' Simulate new data from the Heffernan-Tawn model
#'
#' @param HT a list returned by \code{fit_HT}
#' @param nsim number of simulations
#' @param pth scalar for the quantile probability of the
#' conditioning variable above which to simulate
#' @return a matrix of simulations for the conditioning variable exceeding \code{pth} quantile
sim_HT <- function(HT, nsim = 1e4L, pth){
   stopifnot(length(pth) == 1L, pth > 0, pth < 1)
   u <- qlaplace(pth)
   stopifnot(u > HT$u)
   Z <- HT$resid
   alpha <- HT$alpha
   beta <- HT$beta
   X0.s <- rexp(nsim) + u
   Z.s <- Z[sample(1:nrow(Z), size = nsim, replace = TRUE), , drop=FALSE]
   X.s <- matrix(0, nrow = nsim, ncol = ncol(Z) + 1L)
   for(j in 1:length(alpha)){
      X.s[, j + I(j>=HT$cond.var)] <- alpha[j]*X0.s + X0.s^beta[j]*Z.s[,j]
   }
    X.s[, HT$cond.var] <- X0.s
   return(X.s)
}

#'  Heffernan-Tawn model with fixed margins
#'
#' @param data matrix or data frame of observations on the standard Gumbel scale
#' @param cond.var conditioning variables
#' @param nsim number of simulation samples from the model
#' @param thresh quantile level (uniform scale)
#' @return a vector of length 2 with the probability estimates
probs_HT <- function(data,
                  thresh,
                  nsim = 1e5L) {
   data <- gumbel2laplace(as.matrix(data))
   yp <- pgumbel(6)
   vp <- pgumbel(7)
   v <- qlaplace(vp)
   y <- qlaplace(yp)
   p1 <- mean(sapply(1:3, function(i){
      HT <- fit_HT(data = data, cond.var = i, thresh = thresh)
      simDat <- sim_HT(HT = HT,
                       nsim = nsim,
                       pth = yp)
      mean(rowSums(simDat > y) == 3L)
   })) * (1-yp)
   p2 <- mean(sapply(1:2, function(i){
      HT <- fit_HT(data = data, cond.var = i, thresh = thresh)
      simDat <- sim_HT(HT = HT,
                       nsim = nsim,
                       pth = vp)
      mean(simDat[,1] > v & simDat[,2] > v & simDat[,3] < 0)
   })) * (1 - vp)
 c(p1 = p1, p2 = p2)
}


#' Obtain conditional radius via rolling windows
#'
quantR <- function(data, thresh){
   x <- as.matrix(data)  ## data are already on the standard Gumbel scale
   x <- qexp(pgumbel(x)) ## transform data to standard exponential scales
   # stopifnot(imp.weight >= 1)
   r <- apply(x, 1, sum)
   w <- x / r
   geometricMVE::QR.3d(r = r, w = w, tau = thresh)
}

#' Fit and simulate from the geometric model of Wadsworth and Campbell
#'
#' @param data \code{n} by 3 data matrix data with standard Gumbel marginals
#' @param gauge.fun gauge function, in total three gauge functions are possible. (1) gauge_rvad_full3d, (2) gauge_gaussian3d (3) gauge_clayton3d
#' @param initial.val initial values for the model parameters.  For example, for p=3 and gauge function (1) c(3,0.5), (2) c(3,rep(0.5,len=3)), (3) c(3,0.5)
#' @param lower.limit lower limit of the parameters, usually needed for gauge function (1), e.g., lower.limit = c(0,0)
#' @param upper.limit upper limit of the parameters, usually needed for gauge function (1), e.g., upper.limit = c(100,1)
#' @param nsim  number of simulated example from the fitted model
#' @param thresh quantile level for the radial variable
#'
#' @return
#' @export
#'
#' @return a vector of length 2 with the probability estimates
probs_geometric <-  function(
      data,
      gauge.fun = geometricMVE::gauge_gaussian3d,
      initial.val,
      lower.limit = NULL,
      upper.limit = NULL,
      nsim = 1e5L,
      thresh,
      qr) {
   x <- as.matrix(data)  ## data are already on the standard Gumbel scale
   x <- qexp(pgumbel(x)) ## transform data to standard exponential scales
   y <- 6; v <- 7
   v <- qexp(pgumbel(v))
   y <- qexp(pgumbel(y))
   m <- qexp(0.5)
   # stopifnot(imp.weight >= 1)
   r <- apply(x, 1, sum)
   w <- x / r
   if(missing(qr)){
    qr <- QR.3d(r = r, w = w, tau = thresh)
   }
   rmax <- max(qr$r.tau.wpts, na.rm = TRUE)
   excind <- r > qr$r0w
   rexc <- r[excind]
   wexc <- w[excind, ]

   r0w <- qr$r0w[excind]

   na.ind <- which(is.na(rexc))

   if (length(na.ind) > 0) {
      rexc <- rexc[-na.ind]
      wexc <- wexc[-na.ind, ]
      r0w <- r0w[-na.ind]
   }
   # Fit using different basic parametric gauge functions
   fit <- geometricMVE::fit.geometricMVE.3d(
         r = rexc,
         w = wexc,
         r0w = r0w,
         gfun = gauge.fun,
         init.val = initial.val,
         lower.limit = lower.limit,
         upper.limit = upper.limit,
         control = list(reltol = 1e-6)
      )

   if(max(r0w) > 2*v){
      stop("Threshold level too high")
   }
   ##### Simulate and plot new data above threshold
   xstar <- sim.3d(
         w = wexc,
         r0w = r0w,
         nsim = nsim,
         par = fit$mle,
         gfun = gauge.fun,
         k = 3*y/max(r0w)
      )
      probexc <- (1 - thresh) * Rexc.prob.k.3d(
         k = 3*y/max(r0w),
         r0w = r0w,
         w = wexc,
         gfun = gauge.fun,
         par = fit$mle,
         add.gauge = FALSE
      )
    p1 <- mean(rowSums(xstar > y) == 3L) * probexc
   xstar <- sim.3d(
      w = wexc,
      r0w = r0w,
      nsim = nsim,
      par = fit$mle,
      gfun = gauge.fun,
      k = 2*v/max(r0w)
   )
   #xstar.g2<-sim.3d(w = wexc,r0w = r0w,nsim = nsim,par = fit$mle, gfun = gauge.fun, k=imp.weight) ## k=2, imply that importance weights are used to simulate from the distribution
      probexc <- (1 - thresh) * Rexc.prob.k.3d(
         k = 2*v/max(r0w),
         r0w = r0w,
         w = wexc,
         gfun = gauge.fun,
         par = fit$mle,
         add.gauge = FALSE
      )
   p2 <- mean(xstar[,1] > v & xstar[,2] > v & xstar[,3] < m) * probexc
   c(p1, p2)
}


#' Estimate the two probabilities using the geometric approach followed by the hidden regular variation and regular variation
#'
#' @param data data matrix of dimension \code{n} by 3 with standard Gumbel margins
#' @param quantiles vector of quantiles at which to estimate \eqn{eta}
#' @param a1 probability level for empirical estimation of \code{p1}
#' @param a2 probability level for empirical estimation of \code{p2}
#' @return a list with elements
#' \itemize{
#' \item \code{quantiles}: the vector of quantile levels
#' \item \code{prob} a vector of length 2 giving the average probability estimate at each quantile
#' \item \code{p1} a vector with estimates of \code{p1} at each quantile
#' \item \code{p1} a vector with estimates of \code{p2} at each quantile
#' }
#' @export
probs_hrv <- function(data, quantiles = seq(0.9, 0.99, by = 0.01), a1, a2, 
                      qlev = c("empirical", "theoretical")) {
   qlev <- match.arg(qlev)
   wsE <- qexp(pgumbel(as.matrix(data))) # put in exponential margins
   ye <- qexp(pgumbel(6))
   # qa1 <- qexp(a1)
   Min <- apply(wsE, 1, min)
   #### p1 based on average of different values of eta
   etas_1 <- t(sapply(quantiles, function(q1) {
      thresh <- quantile(Min, q1)
      n <- sum(Min > thresh)
      eta <- mean(Min[Min > thresh] - thresh)
      se.eta <- eta / sqrt(n)
      pmax(0, pmin(eta + c(-1, 0, 1) * 1.96 * se.eta, 1))
   }))
   qa1 <- ifelse("qlev" == "empirical",
                 quantile(Min, probs = a1),
                 qexp(a1))
   p1 <- mean(Min > qa1) * exp(-(ye - qa1) / etas_1[, 2])

   #### probabilities of p2 using p2.star
   data_p2 <-  wsE[data[, 3] < -log(-log(0.5)),1:2]
   ve <- qexp(exp(-exp(-7)))
   # qa2 <- qexp(a2)
   Min <- apply(data_p2, 1, min)
   #### p1 based on average of different values of eta
   etas_2 <- t(sapply(quantiles, function(q1) {
      thresh <- quantile(Min, q1)
      n <- sum(Min > thresh)
      eta <- mean(Min[Min > thresh] - thresh)
      se.eta <- eta / sqrt(n)
      pmax(0, pmin(eta + c(-1, 0, 1) * 1.96 * se.eta, 1))
   }))
   qa2 <- ifelse("qlev" == "empirical",
                 quantile(Min, probs = a2),
                 qexp(a2))
   p2 <- mean(Min > qa2) * exp(-(ve-qa2) / etas_2[, 2]) * 0.5
   ### multiply 1/2 because of the third component less than the median
     return(list(
      quantiles = quantiles,
      "probs" = c(mean(p1), mean(p2)),
      eta1 = etas_1,
      eta2 = etas_2,
      "p1" = p1,
      "p2" = p2
   ))
}

#' Calculating all the probabilities for different models
#'
#' @param data data matrix of dimension n by 3
#' @param mqu marginal quantile probability used in the conditional extremes model, above which to fit a generalized Pareto distribution
#' @param dqu marginal quantile for the conditioning variable of the conditional extremes model
#' @param nsim number of simulated samples to approximate the two probabilities in case of HT and geometricMEV model
#' @param imp.weight whether to use importance sampling to simulate from the fitted model
#' @param thresh quantile level (uniform scale) for conditioning variables for fixed margins for the geometric model (radial threshold), Heffernan-Tawn conditional extremes (marginal threshold) and hidden regular variation (threshold of minimum)
#' @param a1  marginal quantile for empirical extrapolation of the tail probability for p1
#' @param a2  marginal quantile for empirical extrapolation of the tail probability for p2
#'
#' @return
#' @export
#'
#' @examples
estimate_probs_C3 <-
   function(data = Coputopia,
            nsim = 1e6,
            a1 = 0.96,
            a2 = 0.97,
            thresh) {
      # Heffernan-Tawn conditional extremes model, fixed margins
      HT_fixedmarg <- replicate(n = 1000, probs_HT(data = Coputopia, thresh = thresh, nsim = nsim))
      HT_prob <- rowMeans(HT_fixedmarg)
      HT_sd <- apply(HT_fixedmarg, 1, sd)/sqrt(ncol(HT_fixedmarg))
      # Hidden regular variation, semiparametric extrapolation
      HRV <- probs_hrv(quantiles = thresh,
                       data = Coputopia,
                       a1 = a1,
                       a2 = a2)
      qRr <- quantR(data = data, thresh = thresh)
      geom <- replicate(n = 1000, 
                        expr = { probs_geometric(
                           data = data,
                           gauge.fun = gauge_gaussian3d,
                           initial.val = c(1.7, rep(0.33, length.out = 3)),
                           nsim = nsim,
                           thresh = thresh,
                           qr = qRr
                        )}) # Geometric approach (Wadsworth et.al. (2023+))   
      geom_prob <- rowMeans(geom)
      geom_sd <- apply(geom, 1, sd)/sqrt(1000)
      probs.matrix <-
         rbind(
            # probs_texmex_HT(
            #    data = data,
            #    threshold = dqu,
            #    mqu = mqu,
            #    nsim = nsim),
            HT_prob,
            probs_hrv(
               data = data,
               quantiles = thresh,
               a1 = a1,
               a2 = a2
            )$probs,
            geom_prob
         )
      colnames(probs.matrix) <- c("p1", "p2")
      rownames(probs.matrix) <- c("HT", "HRV", "geometric")
      attr(probs.matrix, "sd") <- list(HT = HT_sd, geom = geom_sd)
      return(probs.matrix)
   }

