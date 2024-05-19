# Map data to standard Laplace scale
qlaplace <- function(x) {
  ifelse(x < 0.5, log(2 * x), -log(2 * (1 - x)))
}


# 'texmex' function to check whether constraints are satisfied
texmex_ConstraintsAreSatisfied <- function (a, b, z, zpos, zneg, v) {
  C1e <- a <= min(1, 1 - b * min(z) * v ^ (b - 1), 1 - v ^ (b - 1) * min(z) + min(zpos) /
                    v) &
    a <= min(1, 1 - b * max(z) * v ^ (b - 1), 1 - v ^ (b - 1) * max(z) + max(zpos) /
               v)
  C1o <- a <= 1 &
    a > 1 - b * min(z) * v ^ (b - 1) &
    a > 1 - b * max(z) * v ^ (b - 1) &
    (1 - 1 / b) * (b * min(z)) ^ (1 / (1 - b)) * (1 - a) ^ (-b / (1 - b)) + min(zpos) > 0 &
    (1 - 1 / b) * (b * max(z)) ^ (1 / (1 - b)) * (1 - a) ^ (-b / (1 - b)) + max(zpos) >  0
  C2e <- -a <= min(1, 1 + b * v ^ (b - 1) * min(z), 1 + v ^ (b - 1) * min(z) - min(zneg) /
                     v) &
    -a <= min(1, 1 + b * v ^ (b - 1) * max(z), 1 + v ^ (b - 1) * max(z) - max(zneg) /
                v)
  C2o <- -a <= 1 &
    -a > 1 + b * v ^ (b - 1) * min(z) &
    -a > 1 + b * v ^ (b - 1) * max(z) &
    (1 - 1 / b) * (-b * min(z)) ^ (1 / (1 - b)) * (1 + a) ^ (-b / (1 - b)) - min(zneg) > 0 &
    (1 - 1 / b) *
    (-b * max(z)) ^ (1 / (1 - b)) * (1 + a) ^ (-b / (1 - b)) - max(zneg) >   0
  if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
    warning("Strayed into impossible area of parameter space")
    C1e <- C1o <- C2e <- C2o <- FALSE
  }
  (C1e | C1o) && (C2e | C2o)
}

# Compute constants for nonlinear inequality constraints
data_constr <- function(par,
                        data,
                        thresh = 0.95,
                        index = 1:ncol(data),
                        ...) {
  zpos <- c()
  zneg <- c()
  z <- c()
  v <- c()
  alpha <- par[1]
  beta <- par[2]
  u <- qexp(2 * thresh - 1)
  for (ii in seq_along(index)) {
    i <- index[ii]
    X <- data[data[, i] > u, -i, drop = FALSE]
    X0 <- data[data[, i] > u, i]
    Z <- (X - X0 * alpha) / (X0 ^ beta)
    v <- max(c(v, X0))
    zpos <- range(c(zpos, min(apply(X, 1, min) - X0), max(apply(X, 1, max) - X0)))
    z <- range(c(z, Z))
    zneg <- range(c(zneg, min(apply(X, 1, min) + X0), max(apply(X, 1, max) + X0)))
  }
  list(
    z = z,
    zpos = zpos,
    zneg = zneg,
    v = v + 0.1
  )
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
    min(1, 1 - b * min(z) * v ^ (b - 1), 1 - v ^ (b - 1) * min(z) + min(zpos) /
          v) - a,
    min(1, 1 - b * max(z) * v ^ (b - 1), 1 - v ^ (b - 1) * max(z) + max(zpos) /
          v) - a
  )
  # In 'texmex', the constraint is violated as soon as we have a FALSE,
  # thanks to lazy evaluation - this means that NAs are ignored if the first arguments
  # are FALSE...
  C1o <- min(
    a - (1 - b * min(z) * v ^ (b - 1)),
    a  - (1 - b * max(z) * v ^ (b - 1)),
    (1 - 1 / b) * (b * min(z)) ^ (1 / (1 - b)) * (1 - a) ^ (-b / (1 - b)) + min(zpos),
    (1 - 1 / b) * (b * max(z)) ^ (1 / (1 - b)) * (1 - a) ^ (-b / (1 - b)) + max(zpos)
  )
  C2e <- min(
    1 + b * v ^ (b - 1) * min(z),
    1 + v ^ (b - 1) * min(z) - min(zneg) / v,
    1 + b * v ^ (b - 1) * max(z),
    1 + v ^ (b - 1) * max(z) - max(zneg) / v
  ) + a
  C2o <- min(
    -a - (1 + b * v ^ (b - 1) * min(z)),
    -a - (1 + b * v ^ (b - 1) * max(z)),
    (1 - 1 / b) * (-b * min(z)) ^ (1 / (1 - b)) * (1 + a) ^ (-b / (1 - b)) - min(zneg),
    (1 - 1 / b) * (-b * max(z)) ^ (1 / (1 - b)) * (1 + a) ^ (-b / (1 - b)) - max(zneg)
  )
  C1 <- suppressWarnings(max(c(C1e, C1o), na.rm = TRUE))
  C2 <- suppressWarnings(max(c(C2e, C2o), na.rm = TRUE))
  min(C1, C2)
}



#' Exchangeable Heffernan-Tawn model pseudo log likelihood
#'
#' @param data matrix or data frame of observations on the Laplace scale
#' @param par vector of parameters, \eqn{\alpha}, \eqn{\beta} and nuisance parameters -- either mean and std. dev. \eqn{\mu} and \eqn{\sigma} (normal) or location, scale and slant \eqn{\eta}, \eqn{\omega} and \eqn{\kappa} (skew-normal)
#' @param index vector of indices over which to condition, by default all
#' @param thresh quantile level (uniform scale)
#' @param type estimating equation, either \code{norm} for normal or \code{skewnorm} for skew-normal distribution
#' @param constrain logical; if \code{TRUE}, apply the constrains from Keef, Papastathopoulos and Tawn (2013)
eht_pll <- function(par,
                    data,
                    index = 1:ncol(data),
                    thresh = 0.95,
                    type = c("norm", "skewnorm"),
                    constrain = FALSE,
                    ...) {
  type <- match.arg(type)
  stopifnot(isTRUE(all(index %in% seq_len(ncol(
    data
  )))))
  stopifnot(length(thresh) == 1L, thresh > 0.5, thresh < 1)
  if (abs(par[1]) > 1 | par[2] > 1) {
    return(-1e16)
  }
  data <- as.matrix(data)
  # Define parameters
  alpha <- par[1]
  beta <- par[2]
  loc <- par[3]
  scale <- exp(par[4])
  # For-loop
  u <- qexp(2 * thresh - 1)
  obj <- rep(0, length(index))
  zpos <- c()
  zneg <- c()
  z <- c()
  v <- c()
  for (ii in seq_along(index)) {
    i <- index[ii]
    X <- data[data[, i] > u, -i, drop = FALSE]
    X0 <- data[data[, i] > u, i]
    Z <- (X - X0 * alpha) / (X0 ^ beta)
    scale_mod <- scale * (X0 ^ beta)
    loc_mod <- alpha * X0 + (X0 ^ beta) * loc
    if (type == "norm") {
      for (j in seq_len(ncol(X))) {
        obj[ii] <- obj[ii] +
          sum(dnorm(
            X[, j],
            mean = loc_mod,
            sd = scale_mod,
            log = TRUE
          ))
      }
    } else if (type == "skewnorm") {
      stopifnot(length(par) >= 5)
      slant = par[5]
      for (j in seq_len(ncol(X))) {
        obj[ii] <- obj[ii] +
          sum(sn::dsn(
            X[, j],
            xi = loc_mod,
            omega = scale_mod,
            alpha = slant,
            log = TRUE
          ))
      }
    }
    zpos <- range(c(zpos, min(apply(X, 1, min) - X0), max(apply(X, 1, max) - X0)))
    z <- range(c(z, Z))
    zneg <- range(c(zneg, min(apply(X, 1, min) + X0), max(apply(X, 1, max) + X0)))
    v <- max(c(v, X0))
  }
  if (constrain) {
    const <- texmex_ConstraintsAreSatisfied(
      a = alpha,
      b = beta,
      z = z,
      zpos = zpos,
      zneg = zneg,
      v = v + 0.1
    )
    if (!isTRUE(const)) {
      obj <- -1e20
    }
  }
  return(obj)
}


eht_constraints <- function(alpha_grid = seq(-1, 1, by = 0.1),
                            beta_grid = seq(-3, 1, by = 0.1),
                            data,
                            index = 1:ncol(data),
                            thresh = 0.95) {
  alpha_grid <- alpha_grid[alpha_grid >= -1 & alpha_grid <= 1]
  stopifnot(length(alpha_grid) >= 1)
  beta_grid <- beta_grid[beta_grid <= 1]
  stopifnot(length(beta_grid) >= 1)
  data <- as.matrix(data)
  constrain_matrix <- matrix(NA,
                             nrow = length(alpha_grid),
                             ncol = length(beta_grid))
  for (i in seq_along(alpha_grid)) {
    for (j in seq_along(beta_grid)) {
      constrain_matrix[i, j] <-
        optim_constr_ht(
          par = c(alpha_grid[i], beta_grid[j]),
          constants = data_constr(
            index = index,
            par = c(alpha_grid[i], beta_grid[j]),
            thresh = thresh,
            data = data
          )
        )
    }
  }
  return(list(
    alpha = alpha_grid,
    beta = beta_grid,
    constrain = constrain_matrix
  ))
}

#' Exchangeable Heffernan-Tawn model profile pseudo log likelihood
#'
#' @inheritParams eht_pll
#' @param alpha_grid grid of values for \eqn{\alpha}
#' @param beta_grid grid of values for \eqn{\beta}
eht_pll_profile <- function(alpha_grid = seq(-1, 1, by = 0.1),
                            beta_grid = seq(-3, 1, by = 0.1),
                            data,
                            index = 1:ncol(data),
                            thresh = 0.95,
                            type = c("norm", "skewnorm")) {
  type <- match.arg(type)
  alpha_grid <- alpha_grid[alpha_grid >= -1 & alpha_grid <= 1]
  stopifnot(length(alpha_grid) >= 1)
  beta_grid <- beta_grid[beta_grid <= 1]
  stopifnot(length(beta_grid) >= 1)
  stopifnot(isTRUE(all(index %in% seq_len(ncol(
    data
  )))))
  stopifnot(length(thresh) == 1L, thresh > 0.5, thresh < 1)
  data <- as.matrix(data)
  u <- qexp(2 * thresh - 1)
  profloglik <- matrix(NA,
                       nrow = length(alpha_grid),
                       ncol = length(beta_grid))
  D <- ncol(data)
  for (i in seq_along(alpha_grid)) {
    for (j in seq_along(beta_grid)) {
      alpha <- alpha_grid[i]
      beta <- beta_grid[j]
      Zs <- c()
      Xmin <- c()
      Xmax <- c()
      X0 <- c()
      for (ii in seq_along(index)) {
        k <- index[ii]
        X <- data[data[, k] > u, -k, drop = FALSE]
        X0s <- data[data[, k] > u, k]
        X0 <- c(X0, X0s)
        for (l in seq_len(ncol(X))) {
          Zs <- c(Zs, (X[, l] - X0s * alpha) / (X0s ^ beta))
        }
        Xmin <- c(Xmin, apply(X, 1, min))
        Xmax <- c(Xmax, apply(X, 1, max))
      }
      zpos <- c(min(Xmin - X0), max(Xmax - X0))
      z <- range(Zs)
      zneg <- c(min(Xmin + X0), max(Xmax + X0))
      constraint <- texmex_ConstraintsAreSatisfied(
        a = alpha,
        b = beta,
        z = z,
        zpos = zpos,
        zneg = zneg,
        v = max(X0) + 0.1
      )
      if (isTRUE(constraint)) {
        if (type == "norm") {
          profloglik[i, j] <- -beta * (D - 1) * sum(log(X0)) +
            sum(dnorm(
              Zs,
              mean = mean(Zs),
              sd = sd(Zs),
              log = TRUE
            ))
        } else{
          fit <- sn::selm(Zs ~ 1)
          profloglik[i, j] <- -beta * (D - 1) * sum(log(X0)) + fit@logL
        }
      }
    }
  }
  return(list(
    alpha = alpha_grid,
    beta = beta_grid,
    pll = profloglik
  ))
}

#' Residuals from exchangeable Heffernan-Tawn model
#'
#' @param alpha scale parameter
#' @param beta power parameter
#' @param data matrix of observations
#' @param thresh vector of thresholds on the uniform scale
#' @param testIndep logical; if \code{TRUE}, compute a test of independence using energy statistics for each conditioning variable in turn
#' @param group integer vector for the group; only observations in group 1 are considered considered to build residuals.
residuals_eht <- function(alpha,
                          beta,
                          data,
                          thresh,
                          testIndep = TRUE,
                          group = rep(1L, ncol(data))) {
  p <- ncol(data)
  stopifnot(length(group) == p, isTRUE(all(group %in% 1:2)))
  ng <- sum(group == 1)
  # Reorder so that observations in group 1 are first
  od <- c(which(group == 1), which(group != 1))
  data <- data[, od]
  thresh <- rep(thresh, ng)
  stopifnot(isTRUE(all(thresh > 0, thresh < 1, is.finite(thresh))))
  qlaplace <- function(x) {
    ifelse(x < 0.5, log(2 * x), -log(2 * (1 - x)))
  }
  u <- qlaplace(thresh)
  alpha <- rep(alpha, p)
  beta <- rep(beta, p)

  # Step 1: Identify exceedances and create vectors of residuals
  # Step 2: Compute maximum of each component (see constraint)
  # Step 3: Check constraint
  sdata <- data[apply(data[, 1:ng], 1, function(x) {
    isTRUE(any(x > u))
  }), ]
  nexc <- rowSums(t(sdata[, 1:ng]) >= u)
  res <- matrix(0, nrow = sum(nexc), ncol = p - 1)
  ind <- 0
  pvalindep <- rep(0, ng)
  for (i in seq_len(ng)) {
    exc <- which(sdata[, i] > u[i])
    Y0 <- sdata[exc, i]
    Z <- (sdata[exc, -i] - alpha[i] * Y0) / (Y0 ^ beta[i])
    if (testIndep) {
      pvalindep[i] <- energy::indep.test(Z, Y0, R = 999)$p.value
    }
    res[(ind + 1):(ind + nexc[i]), ] <- Z
    ind <- ind + nexc[i]
  }
  if (!testIndep) {
    pvalindep <- NULL
  }
  return(list(
    pvalindep = pvalindep,
    ng = ng,
    res = res,
    nexc = nexc
  ))
}


#' Prediction from exchangeable Heffernan-Tawn model
#'
#' @param alpha vector of \eqn{\alpha} parameters, recycled if necessary
#' @param beta vector of \eqn{\beta} parameters, recycled if necessary
#' @param data matrix of data
#' @param thresh vector of uniform thresholds
#' @param region risk region on uniform scale
#' @param B number of Monte Carlo replications
#' @param constraint string indicating whether margins are supposed equiprobable
#' or the relative frequency with which a variable is largest is estimated empirically (multinomial)
#' @param type string; should
predict_eht <- function(alpha,
                        beta,
                        data,
                        thresh,
                        region,
                        B = 1e5,
                        nm = 1e4,
                        constraint = c("none", "equiprob")) {
  p <- ncol(data)
  constraint <- match.arg(constraint)
  thresh <- rep(thresh, p)
  stopifnot(isTRUE(all(thresh > 0, thresh < 1, is.finite(thresh))))
  qlaplace <- function(x) {
    ifelse(x < 0.5, log(2 * x), -log(2 * (1 - x)))
  }
  u <- qlaplace(thresh)
  alpha <- rep(alpha, p)
  beta <- rep(beta, p)
  region <- rep(region, p)
  # Step 1: Identify exceedances and create vectors of residuals
  # Step 2: Compute maximum of each component (see constraint)
  # Step 3: Check constraint
  sdata <- data[apply(data, 1, function(x) {
    isTRUE(any(x > u))
  }), ]
  wmax <- table(apply(sdata, 1, which.max))
  nexc <- rowSums(t(sdata) >= u)
  res <- matrix(0, nrow = sum(nexc), ncol = p - 1)
  ind <- 0
  for (i in seq_len(p)) {
    exc <- which(sdata[, i] > u[i])
    res[(ind + 1):(ind + nexc[i]), ] <- (sdata[exc, -i] - alpha[i] * sdata[exc, i]) /
      (sdata[exc, i] ^ beta[i])
    ind <- ind + nexc[i]
  }
  prob <- rep(0, p)
  for (i in seq_len(p)) {
    nsim <- 0
    while (nsim < B) {
      # Generate exceedances above threshold
      y0 <- rexp(n = nm) + region[i]
      #
      newY <- alpha[i] * y0 +
        res[sample.int(n = nrow(res),
                       size = nm,
                       replace = TRUE), , drop = FALSE] *
        y0 ^ beta[i]
      # Check if the conditioning variable is the largest
      condIsMax <- which(y0 > apply(newY, 1, max))
      # Increment the number of simulations completed
      nsim <- nsim + length(condIsMax)
      # Increment probability by the sum of points
      #  falling in the risk region
      prob[i] <- prob[i] + sum(apply(newY[condIsMax, ], 1, function(x) {
        isTRUE(all(x > region[-i]))
      }))
    }
    # Compute the Monte Carlo average of conditional probability
    # multiplied by marginal probability of exceedance (Y0)
    prob[i] <- prob[i] / nsim * (1 - thresh[i])
  }
  pred <- ifelse(constraint == "equiprob", mean(prob), sum(wmax / sum(wmax) *
                                                             prob))
  return(list(
    pred = pred,
    prob = prob,
    res = res,
    nexc = nexc
  ))
}

qexp_level <- function(Zmin, v) {
  sapply(Zmin, function(z) {
    res <- try(uniroot(
      f = function(y0) {
        z - (v - alpha * y0) / y0 ^ beta
      },
      interval = c(v, 1e4)
    )$root, silent =  TRUE)
    if (inherits(res, "try-error")) {
      return(v) # there is no root,
      # so probability of condition is 1
    } else{
      return(res)
    }
  })
}


logp <- function(y0) {
  -min(y0) - log(length(y0)) +
    log(sum(exp(-y0 + min(y0))))
}

#' Fit multivariate generalized Pareto model
#'
fit_mgp <- function(data,
                    thresh,
                    mthresh,
                    type = c("censored", "uncensored"),
                    riskf = c("max", "min", "sum"),
                    model = c("log", "neglog", "br", "xstud"),
                    confint = c("none", "wald", "lr"),
                    level = 0.95,
                    lambdau = NULL) {
  confint <- match.arg(confint)
  model <- match.arg(model)
  type <- match.arg(type)
  opt <- optimize(
    f = function(par) {
      objective_fun_mgp(
        par = par,
        data = data,
        model = model,
        type = type,
        mthresh = mthresh,
        thresh = thresh,
        lambdau = lambdau
      )
    },
    maximum = TRUE,
    lower = 1e-3,
    upper = 1 - 1e-3
  )
  mle <- opt$maximum
  opt$hessian <- -numDeriv::hessian(
    func = function(par) {
      objective_fun_mgp(
        par = par,
        data = data,
        model = model,
        type = type,
        mthresh = mthresh,
        thresh = thresh,
        lambdau = lambdau
      )
    },
    x = mle
  )
  sd_par <- try(sqrt(1 / opt$hessian), silent = TRUE)
  stopifnot(level < 1, level > 0)
  if (confint == "none") {
    conf <- NULL
  } else {
    if (inherits(sd_par, "try-error")) {
      conf <- NULL
      sd_par <- 0.5 / sqrt(nrow(data))
    } else{
      sd_par <- c(sd_par)
      if (confint == "wald") {
        alpha <- (1 - level) / 2
        conf <- as.numeric(mle + c(-1, 1) * qnorm(1 - alpha / 2) * sd_par)
      }
    }
    if (confint == "lr") {
      psi <- pmin(1 - 1e-3, pmax(1e-3, seq(
        from = -max(0.02, 5 * sd_par),
        to = max(0.02, 5 * sd_par),
        length.out = 31
      ))) + mle
      proflik <- sapply(psi, function(par) {
        objective_fun_mgp(
          par = par,
          data = data,
          model = model,
          type = type,
          mthresh = mthresh,
          thresh = thresh,
          lambdau = lambdau
        )
      })
      out_list <- list(
        mle = mle,
        psi = psi,
        psi.max = mle,
        pll = proflik,
        maxpll = opt$objective
      )
      class(out_list) <- "eprof"
      conf <- as.numeric(confint(out_list, level = level)[2:3])
    }
  }
  return(list(
    mle = as.numeric(mle),
    maxll = as.numeric(opt$objective),
    confint = conf
  ))
}

#' Objective function for multivariate generalized Pareto
#'
#' This function is a wrapper that computes the log likelihood
#' for various parametric models using censored or regular likelihood
#' for threshold exceedances and and returns suitable values if
#'  parameters are outside of the range.
objective_fun_mgp <- function(par,
                              data,
                              thresh,
                              mthresh,
                              riskf = c("max", "min", "sum"),
                              type = c("censored", "uncensored"),
                              model = c("log", "neglog", "br", "xstud"),
                              lambdau = NULL) {
  type <- match.arg(type)
  model <- match.arg(model)
  if (model == "xstud") {
    stop("Currently not supported.")
  }
  sdata <- data[apply(data, 1, max) > thresh, ]
  if (is.null(lambdau)) {
    lambdau <- nrow(sdata) / nrow(data)
  } else{
    stopifnot(length(lambdau) == 1L, lambdau > 0, lambdau <= 1)
  }

  D <- ncol(sdata)




  if (model == "log") {
    stopifnot(length(par) == 1L)
    if (par < 0 | par > 1) {
      return(-1e16)
    }
  } else if (model == "neglog") {
    stopifnot(length(par) == 1L)
    if (par < 0) {
      return(-1e16)
    }
  } else if (model == "br") {
    stopifnot(length(par) == 1L)
    if (par < 0 | par > 1) {
      return(-1e16)
    }
    Lambda <- matrix(data = par,
                     nrow = D,
                     ncol = D)
    diag(Lambda) <- 0
  } else if (model == "xstud") {
    stopifnot(length(par) == 2L)
    if (any(par < 0) | par[2] > 1) {
      return(-1e16)
    }
    df <- par[1]
    Sigma <- matrix(par[2], nrow = D, ncol = D)
    diag(Sigma) <- 1
  }
  # if(model != "xstud"){
  if (type == "censored") {
    mev::clikmgp(
      dat = apply(sdata, 1:2, qexp),
      loc = 0,
      scale = 1,
      shape = 0,
      mthresh = rep(qexp(mthresh), D),
      thresh = qexp(thresh),
      par = switch(
        model,
        log = list(alpha = par),
        neglog = list(alpha = par),
        br = list(Lambda = Lambda),
        xstud = list(df = df, Sigma = Sigma)
      ),
      model = model,
      likt = "mgp",
      lambdau = lambdau
    )
  } else{
    mgp::likmgp(
      dat = apply(sdata, 1:2, qexp),
      loc = 0,
      scale = 1,
      shape = 0,
      thresh = qexp(thresh),
      par = switch(
        model,
        log = list(alpha = par),
        neglog = list(alpha = par),
        br = list(Lambda = Lambda),
        xstud = list(df = df, Sigma = Sigma)
      ),
      model = model,
      likt = "mgp",
      lambdau = lambdau
    )


  }
  # This should be the proportion above
  # } else if(model == "xstud"){
  #   corr <- par[2]
  #   df <- par[1]
  #   # The above gives weird results with the 'mev' implementation
  #   # Try instead this mvPot wrapper
  #   stdata <- apply(data, 1:2,
  #                   qgp,
  #                   scale = 1,
  #                   shape = 1,
  #                   loc = 1)
  #   corrfun <- function(h, par = corr){
  #     stopifnot(length(par) == 1L, par >= 0, par < 1, is.numeric(par), is.finite(par))
  #     if(norm(h, type = "2") == 0){1} else{ corr }
  #   }
  #   h <- matrix(1:D, nrow = D, ncol = 2)
  #   p <- 499L
  #   latticeRule <- genVecQMC(p = p, d = (nrow(h) - 1))
  #   primeP <- latticeRule$primeP
  #   vec <- latticeRule$genVec
  #   obs <- split(stdata, row(stdata))
  #   obs <- obs[sapply(obs, max) > mev::qgp(thresh, 1, 1, 1)]
  #   mvPot::censoredLikelihoodXS(
  #     obs = obs,
  #     loc = h,
  #     corrFun = corrfun,
  #     nu = df,
  #     u = rep(qgp(mthresh, 1, 1, 1), length.out = D),
  #     vec = vec)
  # }
}

#' Construct exchangeable variogram matrix
#'
#' @param par scalar parameter for the variogram entries
#' @param D dimension of the resulting matrix
Lmat_exch <- function(par, D) {
  stopifnot(length(par) == 1L)
  Lambda <- matrix(data = par,
                   nrow = D,
                   ncol = D)
  diag(Lambda) <- 0
  return(Lambda)
}

#' Construct compound symmetry correlation matrix
#'
#' @param par scalar correlation parameter
#' @param D dimension of the resulting matrix
compoundsymmetry <- function(par, D) {
  stopifnot(length(par) == 1L)
  Sigma <- matrix(data = par,
                  nrow = D,
                  ncol = D)
  diag(Sigma) <- 1
  return(Sigma)
}


.weightsBRmin <- function(z, Lambda, ...) {
  ellipsis <- list(...)
  D <- length(z)
  stopifnot(ncol(Lambda) == D | nrow(Lambda) == D)
  weights <- rep(0, D)
  for (j in 1:D) {
    weights[j] <- TruncatedNormal::mvNqmc(
      u = rep(Inf, D - 1),
      n = 1e5,
      l = 2 * Lambda[-j, j] + log(z[-j]) - log(z[j]),
      Sig = 2 * (outer(Lambda[-j, j], Lambda[j, -j], FUN = "+") - Lambda[-j, -j])
    )$prob
  }
  return(weights)
}

#' Min risk functional risk region
#'
#' Probability that a Brown-Resnick random vector with variogram
#' \code{Lambda} in standard Pareto margins exceeds \code{z} in all components
expme_BR_min <- function(z, Lambda) {
  weights <- .weightsBRmin(z = z, Lambda = Lambda)
  sum(weights / z)
}

expme_xstud_min <- function(z, Sigma, df, method = "TruncatedNormal") {
  weights <- .weightsXstud_min(z = z, Sigma = Sigma, df = df)
  sum(weights / z)
}

expme_neglog_min <- function(u, par) {
  stopifnot(length(par) == 1, par > 0)
  sum(u) ^ (-1 / par - 1) * sum(u ^ (1 / par))
}

expme_log_min <- function(u, par) {
  # par in (0,1), beta = 1/par
  d <- length(u)
  stopifnot(d > 1)
  if (isTRUE(all.equal(diff(u), rep(0, d - 1)))) {
    sum(1 / u) * sum(sapply(0:(d - 1), function(k) {
      choose(d - 1, k) * (-1) ^ k * (k + 1) ^ (par - 1)
    }))
  } else{
    xiu <- rep(1, d)
    for (j in seq_len(d)) {
      kj <- (u[-j] / u[j]) ^ (-1 / par)
      xiu[j] <- 1 + sum(sapply(seq_along(kj), function(k) {
        ifelse(k %% 2 == 0, 1, -1) * sum((1 + combn(kj, k, FUN = sum)) ^ (par - 1))
      }))
    }
    return(sum(xiu / u))
  }
}


.weightsXstud_min <- function(z,
                              Sigma,
                              df,
                              method = c("mvtnorm", "mvPot", "TruncatedNormal"),
                              ...) {
  method <- match.arg(method, choices = c("mvtnorm", "mvPot", "TruncatedNormal"))[1]
  ellipsis <- list(...)
  D <- nrow(Sigma)
  stopifnot(nrow(Sigma) == length(z))
  weights <- rep(0, D)
  if (method == "mvtnorm") {
    for (j in 1:D) {
      weights[j] <- mvtnorm::pmvt(
        lower = rep(-Inf, D - 1),
        df = df + 1,
        upper = -exp((log(z[-j]) - log(z[j])) / df) + Sigma[-j, j],
        sigma = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE]) / (df + 1)
      )
    }
  } else if (method == "mvPot") {
    if (!is.null(ellipsis$prime)) {
      prime <- ellipsis$prime
      if (!is.null(ellipsis$genvec)) {
        genVec <- ellipsis$genvec
        stopifnot(length(genVec) == D - 1)
      }
    } else {
      prime <- 499L
      genVec <- mvPot::genVecQMC(p = prime, D - 1)$genVec
    }
    for (j in 1:D) {
      weights[j] <- mvPot::mvTProbQuasiMonteCarlo(
        p = prime,
        upperBound = -exp((log(z[-j]) - log(z[j])) / df) + Sigma[-j, j],
        cov = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*% Sigma[j, -j, drop = FALSE]) / (df + 1),
        nu = df + 1,
        genVec = genVec
      )[1]
    }
  } else if (method == "TruncatedNormal") {
    for (j in 1:D) {
      weights[j] <- TruncatedNormal::mvTqmc(
        l = rep(-Inf, D - 1),
        df = df + 1,
        n = 1e5,
        u = -exp((log(z[-j]) - log(z[j])) / df) + Sigma[-j, j],
        Sig = (Sigma[-j, -j] - Sigma[-j, j, drop = FALSE] %*%
                 Sigma[j, -j, drop = FALSE]) / (df + 1)
      )$prob
    }
  }
  return(weights)
}

#' Empirical extrapolation for exchangeable data
#'
#' Compute eta coefficient for all values in the cluster
#' at a lower threshold and use Ledford and Tawn (1997) to extrapolate
#' We do the extrapolation in exponential scale. For unequal levels,
#' we use the angular tail dependence function of Wadsworth and Tawn (2013)
#'
#' @param qlev quantile level
#' @param edat data on exponential margins
#' @param group integer vector indicating the tail probability from \code{stailprob}
#' @param stailprob tail probability at which to compute the probability
tail_emp_exchange <- function(qlev,
                              edata,
                              group = rep(1, ncol(edata)),
                              stailprob = c(1 / 300, 12 / 300)) {
  stopifnot(length(qlev) == 1L, qlev < 1, qlev > 0.5)
  group <- as.integer(group)
  s_exp <- qexp(stailprob[sort(unique(group))], lower.tail = FALSE)
  if (isTRUE(all.equal(group, rep(1, length(group))))) {
    # Structural variable
    Texp <- apply(edata, 1, min)
    # Pick threshold
    s0_exp <- quantile(Texp, qlev)
    eta <- mean(Texp[Texp > s0_exp]) - s0_exp
    # Compute empirical probability that the structure variable min
    # is above upper bound
    # p_ub <- mean(apply(unif[ , which(hclust_od == cl)], 1, min) > ub)

    return(list(
      qlev = as.numeric(qlev),
      eta = as.numeric(eta),
      prob = as.numeric(exp(-(s_exp[1] - s0_exp) / pmin(eta, 1)) * (1 - qlev))
    ))
  } else{
    stopifnot(isTRUE(all(group %in% 1:2)))
    usize <- as.integer(table(group))
    # Compute size of cluster and number of variables in U1 and U2
    u2size <- usize[2]
    clsize <- sum(usize)
    # Under exchangeability, data are exchangeable
    combos <- combn(m = u2size, x = clsize)
    # Estimate eta for weighted data for all permutations
    eta_cl <- tail_prob <- numeric(ncol(combos))
    for (j in 1:ncol(combos)) {
      weights <- rep(1, clsize)
      weights[combos[, j]] <- s_exp[2] / s_exp[1]
      wedat <- sweep(edata, 2, weights, FUN = "/")
      un <- quantile(apply(wedat, 1, min), qlev)
      Texp <- apply(wedat, 1, min)
      eta_cl[j] <- mean(Texp[Texp > un] - un)
      tail_prob[j] <- exp(-(s1_exp - un) / pmin(1, eta_cl[j]))
    }
    return(list(
      qlev = as.numeric(qlev),
      eta = as.numeric(mean(eta_cl)),
      prob = as.numeric(mean(tail_prob) * (1 - qlev))
    ))
  }
}



#' Logistic point process defined on sub-region
#'
#' @param tdat data matrix on unit Pareto scale
#' @param alpha coefficient of the logistic model
#' @param tu threshold vector on the unit Pareto scale
logistic.min.loglik <- function(dat, alpha, u) {
  if (missing(alpha)) {
    stop("Invalid \"par\"")
  }
  alpha <- alpha[1]
  if (alpha > 1) {
    alpha <- 1 / alpha
  }
  if (alpha < 0) {
    stop("Invalid \"par\" for \"log\" model.")
  }
  Xiu <- expme_log_min(u = rep(u, ncol(dat)), par = alpha)
  lVfunlog <- function(x, alpha) {
    if (is.null(dim(x))) {
      alpha * log(sum(x ^ (-1 / alpha)))
    }
    else {
      alpha * log(rowSums(x ^ (-1 / alpha)))
    }
  }
  lVx <- sum(lVfunlog(x = dat, alpha = alpha))
  lfalfacto1 <- function(x, s) {
    sum(log(abs(seq(x, x - s + 1, by = -1))))
  }
  ldVfunlog <- function(x, alpha, lV) {
    falf <- lfalfacto1(alpha, ncol(x))
    - length(dat) * log(alpha) + nrow(x) * falf -
      (1 / alpha + 1) * sum(log(x)) + (alpha - ncol(x)) * lV / alpha
  }
  intens <- ldVfunlog(x = dat, alpha = alpha, lV = lVx)
  intens - nrow(dat) * log(Xiu)
}

neglog.min.loglik <- function(dat, alpha, u) {
  if (missing(alpha)) {
    stop("Invalid \"par\"")
  }
  alpha <- alpha[1]
  if (alpha < 0) {
    alpha <- -alpha
  }
  Xiu <- expme_neglog_min(u = rep(u, ncol(dat)), par = alpha)
  dlogVneglog <- function(x, alpha) {
    x <- as.matrix(x)
    nrow(x) * (lgamma(1 / alpha + ncol(x)) - lgamma(1 / alpha)) +
      length(x) * log(alpha) + (alpha - 1) * sum(log(x)) -
      (1 / alpha + ncol(x)) * sum(log(rowSums(x ^ alpha)))
  }
  intens <- dlogVneglog(x = dat, alpha = alpha)
  intens - nrow(dat) * log(Xiu)
}

hr.min.loglik <- function(dat, alpha, u) {
  stopifnot(length(alpha) == 1L)
  if (alpha < 0 | alpha > 1) {
    return(-1e16)
  }
  D <- ncol(dat)
  Lambda <- matrix(data = alpha,
                   nrow = D,
                   ncol = D)
  diag(Lambda) <- 0
  intens <- mev::intensBR(tdat = dat, Lambda = Lambda)
  Xiu <- expme_BR_min(z = rep(u, D), Lambda = Lambda)
  intens - nrow(dat) * log(Xiu)
}

fit_mingp <- function(data,
                      thresh,
                      family = c("log", "neglog", "hr"),
                      lb = 0.1,
                      ub = 0.99) {
  family <- match.arg(family)
  data <- as.matrix(data[apply(data, 1, min) > thresh, ])
  tdat <- mev::qgp(data,
                   scale = 1,
                   shape = 1,
                   loc = 1)
  tu <- mev::qgp(thresh,
                 scale = 1,
                 shape = 1,
                 loc = 1)
  # Subset only extremes which satisfy the min constraint
  optimfn <- switch(family,
                    log = logistic.min.loglik,
                    neglog = neglog.min.loglik,
                    hr = hr.min.loglik)
  opt <- optimize(
    f = function(alpha) {
      -optimfn(dat = tdat, alpha = alpha, u = tu)
    },
    interval = c(lb, ub)
  )
  loglik <- opt$objective
  mle <- opt$minimum

  lower <- optimize(
    f = function(alpha) {
      abs(2 * (loglik + optimfn(tdat, alpha = alpha, u = tu))
          + qchisq(0.95, 1))
    },
    interval = c(lb, mle)
  )$minimum
  upper <- optimize(
    f = function(alpha) {
      abs(2 * (loglik + optimfn(tdat, alpha = alpha, u = tu))
          + qchisq(0.95, 1))
    },
    interval = c(mle, ub)
  )$minimum
  return(c(
    mle = mle,
    loglik = -loglik,
    lower = lower,
    upper = upper,
    nexc = nrow(data)
  ))
}
