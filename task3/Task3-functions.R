
# Map data to standard Laplace scale
qlaplace <- function(x){
  ifelse(x < 0.5,
         log(2*x),
         -log(2*(1-x))
  )
}

#' Distribution function of the Gumbel distribution
#'
#' @param q a vector or matrix of quantiles
#' @return a vector or matrix of probabilities
pgumbel <- function(q) {
  exp(-exp(-q))
}


#' Empirical estimation of the tail probability
#'
#' @param data \code{n} by 3 data matrix with uniform margins
#' @param levs probability levels passed on uniform scales
#' @param prob string; the probability to estimate, either \code{p1} or \code{p2}.
#'
#' @return estimated probability
#' @export
emp_C3 <- function(
  data,
  levs = pgumbel(c(6, 7,-log(log(2)))),
  # c(6, 7, -log(log(2))),
  prob = c("p1", "p2")) {
  prob <- match.arg(
   prob,
   choices = c("p1", "p2"),
   several.ok = TRUE)
  Y1 <- data[, 1]
  Y2 <- data[, 2]
  Y3 <- data[, 3]
  res <- c()
  if ("p1" %in% prob) {
   res <- c(res, p1 = mean(apply(data, 1, min) > levs[1]))
  }
  if ("p2" %in% prob) {
   res <- c(res, p2 = mean((Y1 > levs[2]) & (Y2 > levs[2]) & (Y3 < levs[3])))
  }
   return(res)
}

# Adaptive barrier; returns negative value when violated
optim_constr_ht <- function (par, constants, ...) {
  a <- par[1]
  b <- par[2]
  z <- constants$z
  zpos <- constants$zpos
  zneg <- constants$zneg
  v <- constants$v
  C1e <- min(
    min(1, 1 - b * min(z) * v^(b - 1), 1 - v^(b -1) * min(z) + min(zpos)/v) - a,
    min(1, 1 - b * max(z) * v^(b - 1), 1 - v^(b - 1) * max(z) + max(zpos)/v) - a
  )
  # In 'texmex', the constraint is violated as soon as we have a FALSE,
  # thanks to lazy evaluation - this means that NAs are ignored if the first arguments
  # are FALSE...
  C1o <- min(
    a - (1 - b * min(z) * v^(b - 1)),
    a - (1 - b * max(z) * v^(b - 1)),
    (1 - 1/b) * (b * min(z))^(1/(1 - b)) * (1 - a)^(-b/(1 - b)) + min(zpos),
    (1 - 1/b) * (b * max(z))^(1/(1 - b)) * (1 - a)^(-b/(1 - b)) + max(zpos)
  )
  C2e <- min(
    1 + b * v^(b - 1) * min(z),
    1 + v^(b - 1) * min(z) - min(zneg)/v,
    1 + b * v^(b - 1) * max(z),
    1 + v^(b - 1) * max(z) - max(zneg)/v) + a
  C2o <- min(
    -a - (1 + b * v^(b - 1) * min(z)),
    -a - (1 + b * v^(b - 1) * max(z)),
    (1 - 1/b) * (-b * min(z))^(1/(1 - b)) * (1 + a)^(-b/(1 - b)) - min(zneg),
    (1 - 1/b) * (-b * max(z))^(1/(1 - b)) * (1 + a)^(-b/(1 - b)) - max(zneg))
  C1 <- suppressWarnings(max(c(C1e, C1o), na.rm = TRUE))
  C2 <- suppressWarnings(max(c(C2e, C2o), na.rm = TRUE))
  min(C1, C2)
}

# Compute constants for nonlinear inequality constraints
data_constr_HT <- function(par, X, X0, ...){
 X <- as.matrix(X)
 v <- max(X0) + 0.1
 Z <- matrix(nrow = nrow(X), ncol = ncol(X))
 co <- numeric(ncol(X))
 for(j in seq_len(ncol(X))){
  alpha <- par[j]
  beta <- par[j+ncol(X)]
  Z[,j] <- (X[,j] - X0*alpha)/(X0^beta)
  zpos <- range(X[,j] - X0)
  zneg <- range(X[,j] + X0)
  z <- range(Z[,j])
  co[j] <- optim_constr_ht(par = c(alpha, beta), constants = list(z = z, zpos = zpos, zneg = zneg, v = v))
 }
 return(co)
}



#' Fit the conditional extremes model
#'
#' @param data data on standard Laplace margins
#' @param cond.var integer giving the index of the conditioning variable
#' @param init.val.HT initial values passed in the optimization of parameter values
#' @param thresh threshold on the probability scale
#'
#' @return a list with elements
#' \itemize{
#' \item \code{cond.var} index of the conditioning variable
#' \item \code{threshold} threshold on the standard Laplace scale
#' \item \code{alpha} vector of parameters
#' \item \code{beta} vector of parameters
#' \item \code{resid} matrix of residuals from the fit
#' }
fit_HT_fixed_margin <- function(data,
                cond.var,
                thresh,
                init.val.HT = c(rep(0.2, 2), rep(0.1, 2))) {
 u <- qlaplace(thresh)
 data <- qlaplace(as.matrix(data))
 stopifnot(length(cond.var) == 1L, cond.var <= ncol(data), cond.var >= 1L)
 X <- data[data[, cond.var] > u, -cond.var, drop = FALSE]
 X0 <- data[data[, cond.var] > u, cond.var]
 # p <- ncol(data) - 1L
 ### objective function to be minimized
 fht_pll <- function(par, X, X0){
  # Define parameters
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  alpha <- par[1:p]
  beta <- par[1:p + p]
  obj <- rep(0, p)
  logX0 <- log(X0)
  for (j in seq_len(p)) {
   Z <- (X[,j] - alpha[j] * X0) / X0^beta[j]
   obj[j] <- -beta[j]*sum(logX0) +
    sum(dnorm(Z, mean = mean(Z), sd = sqrt(n/(n-1))*sd(Z), log = TRUE))
  }
  if (!all(is.finite(obj))) {
   obj <- -1e10
  }
  return(-sum(obj))
 }
 p <- ncol(X)
 ## optimization to find the parameter estimates
 ## TODO break the optimization down in smaller chunks!
 ## we can optimize two parameters at a time
 alpha <- beta <- rep(0, p)
 for(j in seq_len(ncol(X))){
 opt <- Rsolnp::solnp(
  pars = init.val.HT[c(j, p + j)],
  fun = fht_pll,
  LB = c(-1, -1e3),
  UB = rep(1, 2L),
  ineqfun = data_constr_HT,
  ineqLB = 0,
  ineqUB = Inf,
  control = list(trace = FALSE),
  X = X[,j, drop = FALSE],
  X0 = X0
 )
 if (opt$convergence != 0) {
   warning("Optimization did not converge.")
 }
 alpha[j] <- opt$pars[1]
 beta[j] <- opt$pars[2]
 }

 # Obtain residuals
 Z <- matrix(0, nrow = nrow(X), ncol = ncol(X))
 for (j in seq_len(p)) {
   Z[, j] <- (X[, j] - alpha[j] * X0) / (X0 ^ beta[j])
 }
 return(list(
  cond.var = cond.var,
  u = u,
  alpha = alpha,
  beta = beta,
  est.param = c(alpha, beta),
  resid = Z
 ))
}

#' Quantile function of the Laplace distribution
#'
#' @param p a vector or matrix of probabilities
#' @return a vector or matrix of quantiles
qlaplace <- function(p) {
 -sign(p - 0.5) * log(2 * pmin(p, 1 - p))
}

#' Map Gumbel data to standard Laplace scale
#'
#' @param data a vector or matrix of quantiles in standard Gumbel margins
#' @return a vector or matrix of quantiles in standard Laplace margins
gumbel2laplace <- function(data) {
 qlaplace(pgumbel(data))
}

#' Simulate new data from the Heffernan-Tawn model
#'
#' @param HT a list returned by \code{fit_HT_fixed_margin}
#' @param nsim number of simulations
#' @param pth scalar for the quantile probability of the
#' conditioning variable above which to simulate
#' @return a matrix of simulations for the conditioning variable exceeding \code{pth} quantile
sim_HT <- function(HT, nsim = 1e4L, pth) {

 stopifnot(length(pth) == 1L, pth > 0, pth < 1)
 u <- qlaplace(pth)
 stopifnot(u > HT$u)
 Z <- HT$resid
 alpha <- HT$alpha
 beta <- HT$beta
 X0.s <- rexp(nsim) + u
 Z.s <- Z[sample(1:nrow(Z), size = nsim, replace = TRUE), , drop = FALSE]
 X.s <- matrix(0, nrow = nsim, ncol = ncol(Z) + 1L)
 for (j in 1:length(alpha)) {
  X.s[, j + I(j >= HT$cond.var)] <-
   alpha[j] * X0.s + X0.s ^ beta[j] * Z.s[, j]
 }
 X.s[, HT$cond.var] <- X0.s
 return(X.s)
}

#' Heffernan-Tawn model with fixed margins
#'
#' @param data matrix or data frame of observations on uniform margins
#' @param nsim number of simulation samples from the model
#' @param init.val.HT initial values for the parameters
#' @param levs probability levels in standard uniform scale
#' @param repl number of replications while simulation the data from the model
#' @param thresh quantile level (uniform scale)
#'
#' @return a vector of length 2 with the probability estimates
probs_HT <- function(data,
           thresh,
           nsim = 1e5L,
           init.val.HT,
           levs = pgumbel(c(6, 7,-log(log(2)))),
           repl = 1L) {
 p1 <- mean(sapply(1:3, function(i) {
  HT <- fit_HT_fixed_margin(
    data = data,
    cond.var = i,
    thresh = thresh,
    init.val.HT = init.val.HT
   )
  cal_prob <- function(fit.HT) {
   simDat <- sim_HT(HT = fit.HT,
            nsim = nsim,
            pth = levs[1])
   return(mean(apply(simDat,1, min) > qlaplace(levs[1])))
  }
  mean(replicate(repl, cal_prob(fit.HT = HT)))
 })) * (1 - levs[1])


 p2 <- mean(sapply(1:2, function(i) {
  HT <-
   fit_HT_fixed_margin(
    data = data,
    cond.var = i,
    thresh = thresh,
    init.val.HT = init.val.HT
   )
  cal_prob <- function(fit.HT) {
   simDat <- sim_HT(HT = fit.HT,
            nsim = nsim,
            pth = levs[2])
   return(simDat[, 1] > qlaplace(levs[2]) &
          simDat[, 2] > qlaplace(levs[2]) &
          simDat[, 3] < qlaplace(levs[3]))
  }
  mean(replicate(repl, cal_prob(fit.HT = HT)))
 })) * (1 - levs[2])

 # m is median, so zero for Laplace distribution
 sdata <- data[data[,3] < levs[3], 1:2]
 p2.star <- mean(sapply(1:2, function(i) {
  HT <- fit_HT_fixed_margin(
    data = sdata,
    cond.var = i,
    thresh = thresh,
    init.val.HT = init.val.HT[c(1, 3)]
   )
  cal_prob <- function(fit.HT) {
   simDat <- sim_HT(HT = fit.HT,
            nsim = nsim,
            pth = levs[2])
   return(mean(apply(simDat, 1, min) > qlaplace(levs[2])))
  }
  mean(replicate(repl, cal_prob(fit.HT = HT)))})) * (1 - levs[2]) * levs[3]
  c(p1, p2, p2.star)
}



#' Fit and simulate from the geometric model of Wadsworth and Campbell
#'
#' @param data \code{n} by 3 data matrix data with standard Gumbel marginals
#' @param gauge.fun gauge function, in total three gauge functions are possible. (1) gauge_rvad_all3d, (2) gauge_gaussian3d (3) gauge_clayton3d
#' @add.gauge If TRUE, then we additively mix (with numerical rescaling) the gauge functions specified by gauge1 and gauge2
#' @param gauge1 ### first gauge function in the mixing gauge function. It could be either of the three enlisted in @param gauge.fun
#' @param gauge2 ### second gauge function in the mixing gauge function. It could be either of the three enlisted in @param gauge.fun
#' @param initial.val initial values for the model parameters. For example, for p=3 and gauge function (1) c(3,0.5), (2) c(3,rep(0.5,len=3)), (3) c(3,0.5)
#' @param lower.limit lower limit of the parameters, usually needed for gauge function (1), e.g., lower.limit = c(0,0)
#' @param upper.limit upper limit of the parameters, usually needed for gauge function (1), e.g., upper.limit = c(100,1)
#' @param nsim number of simulated example from the fitted model
#' @param thresh quantile level for the radial variable
#' @param qr
#'
#' @return
#' @export
#'
#' @return a vector of length 2 with the probability estimates
fit_geometric.3d <- function(data,
               gauge.fun = geometricMVE::gauge_rvad_all3d,
               add.gauge = FALSE,
               gauge1 = NULL,
               gauge2 = NULL,
               initial.val,
               lower.limit = NULL,
               upper.limit = NULL,
               nsim = 1e5L,
               thresh,
               par1ind = 2:length(initial.val),
               par2ind = NULL,
               diagnostic = FALSE) {
 x <- qexp(as.matrix(data)) ## transform data to standard exponential
 r <- apply(x, 1, sum)
 w <- x / r
 qr <- QR.3d(r = r, w = w, tau = thresh)
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

 if (!add.gauge) {
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
   add.gauge = TRUE,
   gauge1 = gauge1,
   gauge2 = gauge2,
   init.val = initial.val,
   lower.limit = lower.limit,
   upper.limit = upper.limit,
   par1ind = par1ind,
   par2ind = par2ind,
   control = list(reltol = 1e-6)
  )
 }

 AIC <- 2 * fit$nllh + 2 * length(fit$mle)
 if(isTRUE(diagnostic)){
  if(!add.gauge){
  geometricMVE::qqdiag.3d(
   r = rexc,
   w = wexc,
   r0w = r0w,
   par = fit$mle,
   gfun = gauge.fun,
   add.gauge = FALSE)

  } else{
   geometricMVE::qqdiag.3d(
    r = rexc,
    w = wexc,
    r0w = r0w,
    par = fit$mle,
    gauge1 = gauge1,
    gauge2 = gauge2,
    par1ind = par1ind,
    par2ind = par2ind,
    add.gauge = TRUE)
  }
 }
 return(list(
  AIC = AIC,
  r0w = r0w,
  mle = fit$mle,
  wexc = wexc,
  gauge.fun = gauge.fun
 ))
}


# Fit geometric extremes model to bivariate data
fit_geometric.2d <- function(data,
               gauge.fun = geometricMVE::gauge_rvad,
               add.gauge = FALSE,
               gauge1 = NULL,
               gauge2 = NULL,
               initial.val,
               lower.limit = NULL,
               upper.limit = NULL,
               nsim = 1e5L,
               thresh,
               par1ind = 2:length(initial.val),
               par2ind = NULL,
               diagnostic = FALSE) {
 ## transform data to standard exponential  
 x <- qexp(as.matrix(data))
 ## compute pseudo-polar coordinates
 r <- rowSums(x)
 w <- x[, -ncol(data), drop = FALSE] / r
 ## obtain thresholds using moving windows
 qr <- QR.2d(
  r = r,
  w = w,
  tau = thresh,
  method = "empirical"
 )
 ## compute exceedances
 rmax <- max(qr$r.tau.wpts, na.rm = TRUE)
 excind <- r > qr$r0w
 rexc <- r[excind]
 wexc <- w[excind]
 r0w <- qr$r0w[excind]

 na.ind <- which(is.na(rexc))
 if (length(na.ind) > 0) {
  rexc <- rexc[-na.ind]
  wexc <- wexc[-na.ind, ]
  r0w <- r0w[-na.ind]
 }

 if (!add.gauge) {
  # Fit using different basic parametric gauge functions
  fit <- geometricMVE::fit.geometricMVE.2d(
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
  fit <- geometricMVE::fit.geometricMVE.2d(
   r = rexc,
   w = wexc,
   r0w = r0w,
   add.gauge = TRUE,
   gauge1 = gauge1,
   gauge2 = gauge2,
   init.val = initial.val,
   lower.limit = lower.limit,
   upper.limit = upper.limit,
   par1ind = par1ind,
   par2ind = par2ind,
   control = list(reltol = 1e-6)
  )
 }

 AIC <- 2 * fit$nllh + 2 * length(fit$mle)
 if(isTRUE(diagnostic)){
  if(!add.gauge){
   geometricMVE::qqdiag.2d(
    r = rexc,
    w = wexc,
    r0w = r0w,
    par = fit$mle,
    gfun = gauge.fun,
    add.gauge = FALSE)

  } else{
   geometricMVE::qqdiag.2d(
    r = rexc,
    w = wexc,
    r0w = r0w,
    par = fit$mle,
    gauge1 = gauge1,
    gauge2 = gauge2,
    par1ind = par1ind,
    par2ind = par2ind,
    add.gauge = TRUE)
  }
 }
 return(list(
  AIC = AIC,
  r0w = r0w,
  mle = fit$mle,
  wexc = wexc,
  gauge.fun = gauge.fun
 ))
}



#' Fit and simulate from the geometric model of Wadsworth and Campbell
#'
#' @param data \code{n} by 3 data matrix data with standard Gumbel marginals
#' @param gauge.fun gauge function, in total three gauge functions are possible. (1) gauge_rvad_full3d, (2) gauge_gaussian3d (3) gauge_clayton3d
#' @param initial.val initial values for the model parameters. For example, for p=3 and gauge function (1) c(3,0.5), (2) c(3,rep(0.5,len=3)), (3) c(3,0.5)
#' @param lower.limit lower limit of the parameters, usually needed for gauge function (1), e.g., lower.limit = c(0,0)
#' @param upper.limit upper limit of the parameters, usually needed for gauge function (1), e.g., upper.limit = c(100,1)
#' @param nsim number of simulated example from the fitted model
#' @param thresh quantile level for the radial variable
#' @export
#'
#' @return a vector of length 2 with the probability estimates
probs_geometric.3d <- function(data,
                nsim = 1e5,
                gfun = geometricMVE::gauge_rvad_all3d,
                k0.1 = 1,
                k0.2 = 1,
                initial.val = c(3, rep(0.33, length.out = 3), 0.9),
                levs = pgumbel(c(6, 7,-log(log(2)))),
                lower.limit = c(0, rep(0, 4)),
                upper.limit = c(5, rep(1, 4)),
                thresh = 0.95,
                add.gauge = FALSE,
                gauge1 = NULL,
                gauge2 = NULL,
                par1ind,
                par2ind,
                repl) {
 x <- qexp(as.matrix(data))

 fit <- fit_geometric.3d(
  data = data,
  gauge.fun = gfun,
  add.gauge = add.gauge,
  gauge1 = gauge1,
  gauge2 = gauge2,
  initial.val = initial.val,
  lower.limit = lower.limit,
  upper.limit = upper.limit,
  nsim = 1e5L,
  thresh = thresh,
  par1ind = par1ind,
  par2ind = par2ind
 )


 cal_prob <- function(fit) {
  ### Estimating probability p1
  if (!add.gauge) {
   xstar <- sim.3d(
    w = fit$wexc,
    r0w = fit$r0w,
    nsim = nsim,
    par = fit$mle,
    gfun = gfun,
    k = k0.1,
    add.gauge = FALSE
   )
   probexc <- (1 - thresh) * Rexc.prob.k.3d(
    k = k0.1,
    r0w = fit$r0w,
    w = fit$wexc,
    gfun = gfun,
    par = fit$mle,
    add.gauge = FALSE
   )
  } else{
   xstar <- sim.3d(
    w = fit$wexc,
    r0w = fit$r0w,
    nsim = nsim,
    par = fit$mle,
    gfun = gfun,
    k = k0.1,
    add.gauge = TRUE,
    gauge1 = gauge1,
    gauge2 = gauge2
   )
   probexc <- (1 - thresh) * Rexc.prob.k.3d(
    k = k0.1,
    r0w = fit$r0w,
    w = fit$wexc,
    gfun = gfun,
    par = fit$mle,
    add.gauge = TRUE,
    gauge1 = gauge1,
    gauge2 = gauge2
   )
  }
   p1 <- mean(rowSums(xstar > qexp(levs[1])) == 3L) * probexc
  ### Estimating probability p2
  if (!add.gauge) {
   xstar <- sim.3d(
    w = fit$wexc,
    r0w = fit$r0w,
    nsim = nsim,
    par = fit$mle,
    gfun = gfun,
    k = k0.1,
    add.gauge = FALSE
   )
   probexc <- (1 - thresh) * Rexc.prob.k.3d(
    k = k0.1,
    r0w = fit$r0w,
    w = fit$wexc,
    gfun = gfun,
    par = fit$mle,
    add.gauge = FALSE
   )
  } else{
   xstar <- sim.3d(
    w = fit$wexc,
    r0w = fit$r0w,
    nsim = nsim,
    par = fit$mle,
    gfun = gfun,
    k = k0.1,
    add.gauge = TRUE,
    gauge1 = gauge1,
    gauge2 = gauge2
   )
   probexc <- (1 - thresh) * Rexc.prob.k.3d(
    k = k0.1,
    r0w = fit$r0w,
    w = fit$wexc,
    gfun = gfun,
    par = fit$mle,
    add.gauge = TRUE,
    gauge1 = gauge1,
    gauge2 = gauge2
   )
  }
   p2 <- mean(xstar[, 1] > qexp(levs[2]) & xstar[, 2] > qexp(levs[1]) & xstar[, 3] < qexp(levs[3])) * probexc
  c(p1, p2)
 }

 geom <- replicate(n = repl, cal_prob(fit = fit))
 probs <- rowMeans(geom)
 return(probs)
}




#' Fit and simulate from the geometric model of Wadsworth and Campbell
#'
#' @param data \code{n} by 3 data matrix data with standard Gumbel marginals
#' @param gauge.fun gauge function, in total three gauge functions are possible. (1) gauge_rvad_full3d, (2) gauge_gaussian3d (3) gauge_clayton3d
#' @param initial.val initial values for the model parameters. For example, for p=3 and gauge function (1) c(3,0.5), (2) c(3,rep(0.5,len=3)), (3) c(3,0.5)
#' @param lower.limit lower limit of the parameters, usually needed for gauge function (1), e.g., lower.limit = c(0,0)
#' @param upper.limit upper limit of the parameters, usually needed for gauge function (1), e.g., upper.limit = c(100,1)
#' @param nsim number of simulated example from the fitted model
#' @param thresh quantile level for the radial variable
#' @export
#'
#' @return a vector of length 2 with the probability estimates
probs_geometric.2d <- function(data,
                nsim = 1e5,
                gfun = geometricMVE::gauge_rvad,
                k0.1 = 1,
                k0.2 = 1,
                initial.val = c(3, 0.9),
                levs = pgumbel(c(6, 7,-log(log(2)))),
                lower.limit = c(0, rep(0, 4)),
                upper.limit = c(5, rep(1, 4)),
                thresh = 0.95,
                add.gauge = FALSE,
                gauge1 = NULL,
                gauge2 = NULL,
                par1ind,
                par2ind,
                repl) {

 fit <- fit_geometric.2d(
  data = data[data[,3] < levs[3],-3],
  gauge.fun = gfun,
  add.gauge = add.gauge,
  gauge1 = gauge1,
  gauge2 = gauge2,
  initial.val = initial.val,
  lower.limit = lower.limit,
  upper.limit = upper.limit,
  nsim = 1e5L,
  thresh = thresh,
  par1ind = par1ind,
  par2ind = par2ind
 )


 cal_prob <- function(fit) {
  ### Estimating probability p2
  if (!add.gauge) {
   xstar <- sim.2d(
    w = fit$wexc,
    r0w = fit$r0w,
    nsim = nsim,
    par = fit$mle,
    gfun = gfun,
    k = k0.2,
    add.gauge = FALSE
   )

   probexc <- (1 - thresh) * Rexc.prob.k.2d(
    k = k0.2,
    r0w = fit$r0w,
    w = fit$wexc,
    gfun = gfun,
    par = fit$mle,
    add.gauge = FALSE
   )
  } else{
   xstar <- sim.2d(
    w = fit$wexc,
    r0w = fit$r0w,
    nsim = nsim,
    par = fit$mle,
    gfun = gfun,
    k = k0.2,
    add.gauge = TRUE,
    gauge1 = gauge1,
    gauge2 = gauge2
   )

   probexc <- (1 - thresh) * Rexc.prob.k.2d(
    k = k0.2,
    r0w = fit$r0w,
    w = fit$wexc,
    gfun = gauge.fun,
    par = fit$mle,
    add.gauge = TRUE,
    gauge1 = gauge1,
    gauge2 = gauge2
   )
  }
   mean(apply(xstar, 1, min) > qexp(levs[2])) * probexc * levs[3]
 }

 geom <- replicate(n = repl, cal_prob(fit = fit))
 probs <- mean(geom)
 probs
}


#' Estimate the two probabilities using the geometric approach followed by the hidden regular variation and regular variation
#'
#' @param data data matrix of dimension \code{n} by 3 with standard Gumbel margins
#' @param quantiles vector of quantiles at which to estimate \eqn{eta}
#' @param a1 probability level for empirical estimation of \code{p1}
#' @param qlev
#' @param levs probability at standard uniform scales for calculation of probabilities
#' @param a2 probability level for empirical estimation of \code{p2}
#'
#' @return a list with elements
#' \itemize{
#' \item \code{quantiles}: the vector of quantile levels
#' \item \code{prob} a vector of length 2 giving the average probability estimate at each quantile
#' \item \code{p1} a vector with estimates of \code{p1} at each quantile
#' \item \code{p1} a vector with estimates of \code{p2} at each quantile
#' }
#' @export
probs_hrv <-
 function(data,
      quantiles = seq(0.9, 0.99, by = 0.01),
      a1,
      a2,
      qlev = c("empirical", "theoretical"),
      levs = pgumbel(c(6, 7,-log(log(2))))) {
  stopifnot(a1 < min(levs[1:2]), a2 < min(levs[1:2]))
  qlev <- match.arg(qlev)
  edata <- qexp(as.matrix(data))
  ye <-  qexp(levs[1])
  # qa1 <- qexp(a1)
  Min <- apply(edata, 1, min)
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
  stopifnot(qa1 < ye)
  p1 <- mean(Min > qa1) * exp(-(ye - qa1) / etas_1[, 2])

  #### probabilities of p2 using p2.star
  edata_p2 <- edata[data[, 3] < levs[3], 1:2]
  ve <- qexp(levs[2])
  # qa2 <- qexp(a2)
  Min <- apply(edata_p2, 1, min)
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
  stopifnot(qa2 < ve)
  p2 <- mean(Min > qa2) * exp(-(ve - qa2) / etas_2[, 2]) * levs[3]
  ### multiply by 1/2 because of the third component less than the median
  return(list(
   quantiles = quantiles,
   "probs" = c(mean(p1), NA, mean(p2)),
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
#' @param a1 marginal quantile for empirical extrapolation of the tail probability for p1
#' @param a2 marginal quantile for empirical extrapolation of the tail probability for p2
#'
#' @return
#' @export
#'
#' @examples
estimate_probs_C3 <- function(
   data,
   nsim = 1e6,
   a1 = 0.96,
   a2 = 0.97,
   thresh,
   init.val.geom.3d = c(3, rep(0.33, 3L), 0.9),
   init.val.geom.2d = c(3, 0.9),
   init.val.HT = c(rep(0.2, 2), rep(0.1, 2)),
   levs = pgumbel(c(6, 7, -log(log(2)))),
   repl = 1L) {
  # Heffernan-Tawn conditional extremes model, fixed margins
  HT_prob <- probs_HT(
   data = data,
   thresh = thresh,
   nsim = 1e5L,
   init.val.HT = init.val.HT,
   levs = levs,
   repl = repl
  )
  # Hidden regular variation, semiparametric extrapolation
  HRV_prob <- probs_hrv(
   quantiles = thresh,
   data = data,
   a1 = a1,
   a2 = a2,
   qlev = "empirical",
   levs = levs
  )$probs
  ### # Geometric approach (Wadsworth and Campbell (2024+))
  geom_prob <- c(
   probs_geometric.3d(
    data = data,
    nsim = nsim,
    gfun = geometricMVE::gauge_rvad_all3d,
    k0.1 = 2,
    k0.2 = 2,
    initial.val = init.val.geom.3d,
    levs = levs,
    lower.limit = c(0, rep(0, 4)),
    upper.limit = c(5, rep(1, 4)),
    thres = thresh,
    repl = repl
   ),
   probs_geometric.2d(
    data = data,
    nsim = nsim,
    gfun = geometricMVE::gauge_rvad,
    k0.1 = 1.5,
    k0.2 = 1.5,
    initial.val = init.val.geom.2d,
    levs = levs,
    lower.limit = c(0, 0),
    upper.limit = c(5, 1),
    thres = thresh,
    repl = repl
   )
  )
  probs.matrix <-
   rbind(HT_prob,
      HRV_prob,
      geom_prob)
  colnames(probs.matrix) <- c("p1", "p2", "p2.star")
  rownames(probs.matrix) <- c("HT", "HRV", "geometric")
  return(probs.matrix)
 }
