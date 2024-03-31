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
#' @add.gauge  If TRUE, then we additively mix (with numerical rescaling) the gauge functions specified by gauge1 and gauge2
#' @param gauge1 ### first gauge function in the mixing gauge function. It could be either of the three enlisted in @param gauge.fun 
#' @param gauge2 ### second gauge function in the mixing gauge function. It could be either of the three enlisted in @param gauge.fun 
#' @param initial.val initial values for the model parameters.  For example, for p=3 and gauge function (1) c(3,0.5), (2) c(3,rep(0.5,len=3)), (3) c(3,0.5)
#' @param lower.limit lower limit of the parameters, usually needed for gauge function (1), e.g., lower.limit = c(0,0)
#' @param upper.limit upper limit of the parameters, usually needed for gauge function (1), e.g., upper.limit = c(100,1)
#' @param nsim  number of simulated example from the fitted model
#' @param thresh quantile level for the radial variable
#' @param qr 
#'
#' @return
#' @export
#'
#' @return a vector of length 2 with the probability estimates
probs_geometric <-  function(
    data,
    gauge.fun = geometricMVE::gauge_gaussian3d, 
    gauge1, 
    gauge2, 
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
  if(!add.gauge){
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
  } else{
    fit <- geometricMVE::fit.geometricMVE.3d(
      r = rexc,
      w = wexc,
      r0w = r0w,
      gauge1 = gauge1,
      gauge2 = gauge2,
      init.val = initial.val,
      lower.limit = lower.limit,
      upper.limit = upper.limit,
      control = list(reltol = 1e-6)
    ) 
  }
  
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
  
  probs<- c(p1, p2)
  AIC<- fit$nllh + 2*length(fit1$mle)
  
  return(list("probs"=probs,  "AIC"=AIC))
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

