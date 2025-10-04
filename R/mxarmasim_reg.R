#' Simulation from the Maxwell ARMA Model with Regressors (MxARMAX)
#'
#' Generates a univariate time series from a Maxwell ARMA-type model that can also
#' incorporate exogenous covariates. Depending on the parameters provided, the
#' function allows simulation of MxARX, MxMAX, or ARMAX-type processes with
#' Maxwell-distributed innovations.
#'
#' @param n Integer. Number of observations to simulate.
#' @param alpha Numeric. Intercept term in the linear predictor.
#' @param phi Numeric vector. Autoregressive coefficients. If \code{NULL},
#'   no AR component is included.
#' @param theta Numeric vector. Moving average coefficients. If \code{NULL},
#'   no MA component is included.
#' @param beta Numeric vector. Regression coefficients associated with the
#'   covariate matrix \code{X}. Must have length equal to the number of columns of \code{X}.
#' @param X Matrix or vector of covariates (regressors). Should have at least
#'   \code{n + 2*max(p,q)} rows, where \code{p = length(phi)} and \code{q = length(theta)}.
#'
#' @details
#' Let \eqn{Y_t} be the simulated process with conditional intensity \eqn{\mu_t}
#' linked to a linear predictor \eqn{\eta_t} through the logarithmic link:
#' \deqn{\mu_t = \exp(\eta_t).}
#' Then,
#' \deqn{\eta_t = \alpha + \mathbf{X}_t^\top\beta +
#'       \sum_{i=1}^{p} \phi_i\left\{ \log(Y_{t-i}) - \mathbf{X}_{t-i}^\top\beta \right\} +
#'       \sum_{j=1}^{q} \theta_j e_{t-j},}
#' where \eqn{e_t = \log(Y_t) - \eta_t} is the one-step-ahead prediction error.
#' The innovations \eqn{Y_t} are generated from the Maxwell distribution via
#' \code{rmax(mu)}.
#'
#' To mimic ARX, MAX, or ARMAX structures, specify only \code{phi}, only
#' \code{theta}, or both, respectively.
#'
#' @return A numeric vector with \code{n} simulated observations.
#'
#' @author Arthur Hintz
#'
#' @examples
#' X  <- stats::runif(1004, 0, 1)
#' y1 <- mxarmareg.sim(n = 1000, alpha = 0.6, phi = c(0.6, 0.1),
#'                     beta = 0.7, X = X)          # ARX
#' y2 <- mxarmareg.sim(n = 1000, alpha = 0.6, theta = 0.3,
#'                     beta = 0.7, X = X)          # MAX
#' y3 <- mxarmareg.sim(n = 1000, alpha = 0.6, phi = c(0.6,0.1),
#'                     theta = 0.2, beta = 0.7, X = X)  # ARMAX
#'
#' @importFrom stats make.link
#'
#' @export
mxarmareg.sim <- function(n, phi = NULL, theta = NULL, alpha = 0.0, beta = NULL, X = NULL) {

    if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive number.")
  }
  if (is.null(phi) && is.null(theta)) {
    stop("At least phi or theta must be specified.")
  }

  #if (is.na(X) || nrow(X) != m + n) {
  # stop("A matriz X deve ter n + m linhas.")}


  if (!is.null(X)) {
    X <- as.matrix(X)}

  is_phi <- !is.null(phi)
  is_theta <- !is.null(theta)

  link <- make.link("log")
  linkfun <- link$linkfun
  linkinv <- link$linkinv


  ## ARX ___________________________________________
  if (is_phi && !is_theta) {

    ar <- 1:length(phi)
    p <- max(ar)
    m <- 2*p

    ynew <- rep(alpha,(n+m))
    mu <- linkinv(ynew)

    eta <- y <- rep(NA, n+m)

    for(i in (m+1):(n+m)){
      eta[i] <- alpha + X[(i),]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-(X[i-ar,]%*%as.matrix(beta))))
      mu[i]   <- linkinv(eta[i])
      y[i]    <- MxARMA::rmax(mu = mu[i])
      ynew[i] <- linkfun(y[i])
    }

    return(y[(m+1):(n+m)])
  }


  ## MAX ________________________________________________
  if (!is_phi && is_theta){

    ma <- 1:length(theta)
    q <- max(ma)
    m <- 2*q

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    eta <- y <- error <- rep(0,n+m)

    for(i in (m+1):(n+m)){
      eta[i] <- alpha + X[(i),]%*%as.matrix(beta) + (theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- MxARMA::rmax(mu = mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]

    }

    return(y[(m+1):(n+m)])
  }

  #ARMAX _______________________________________________
  if(is_phi && is_theta){

    ar <- 1:length(phi)
    ma <- 1:length(theta)
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    error <- rep(0,n+m)
    eta <- y <- rep(NA, n+m)


    for(i in (m+1):(n+m)){
      eta[i]  <- alpha + X[(i),] %*% as.matrix(beta) +  (phi %*% (ynew[i-ar] - X[(i-ar),]%*%as.matrix(beta)))  + theta%*%error[i-ma]
      mu[i]   <- linkinv(eta[i])
      y[i]    <- MxARMA::rmax(mu = mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
    }

    return(y[(m+1):(n+m)])
  }
}
