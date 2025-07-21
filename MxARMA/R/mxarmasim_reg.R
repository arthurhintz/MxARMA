#' Simulation of the Maxwell Distribution with regressors
#'
#' This function generates random values from the MxARMA model with regressors.
#'
#' @param n Number of random values to generate.
#' @param alpha A value for the intercept parameter.
#' @param phi A value or a vector of values for the autoregressive parameter phi.
#' @param theta A value or a vector of values for the moving average parameter theta.
#' @param beta A value or a vector of values for the regressor coefficient.
#' @param X A matrix of covariates or regressors used in the model. It should have the appropriate dimensions for the computations.
#'
#'
#' @details The function simulates values from the MxARMA model, including the possibility of generating values for an autoregressive with exogenous variables (ARX), moving average (MAX), or combined autoregressive moving average with exogenous variables (ARMAX) model.
#'
#' @examples
#'
#'  X = stats::runif(1004, 0, 1)
#'  y <- mxarmareg.sim(n = 1000, alpha = 0.6, phi = c(0.6, 0.1), theta = 0.3, beta = 0.7, X = X)
#'  y
#'
#' @importFrom stats make.link
#'
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
      y[i]    <- MxARMA::rmax(mu[i])
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
      y[i]    <- MxARMA::rmax(mu[i])
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
      y[i]    <- MxARMA::rmax(mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
    }

    return(y[(m+1):(n+m)])
  }
}
