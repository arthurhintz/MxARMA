#' Simulation of the Maxwell Distribution
#'
#' This function generates random values from the MxARMA model.
#'
#' @param n Number of random values to generate.
#' @param alpha A value for the intercept parameter.
#' @param phi A value or a vector of values for the regressor parameter phi.
#' @param theta A value or a vector of values for the moving average parameter theta.
#'
#' @details The function simulates values from the MxARMA model, which includes
#'  the possibility of generating values for an autoregressive (AR), moving average (MA),
#'   or combined autoregressive moving average (ARMA) model.
#'
#' @examples
#'  y <- mxarma.sim(n = 1000, alpha = 0.6, phi = c(0.6, 0.1), theta = 0.3)
#'  y
#'
#' @import stats
#'
#' @export
mxarma.sim <- function(n, alpha = 0.0, phi = NULL, theta = NULL) {

  if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive number.")
  }
  if (!is.numeric(phi) && !is.numeric(theta)) {
    stop("phi or theta must be numeric vectors.")
  }
  if (any(is.null(phi)) && any(is.null(theta))) {
    stop("At least phi or theta must be specified.")
  }

  isphi   <- !is.null(phi)
  istheta <- !is.null(theta)

  link <- stats::make.link("log")
  linkfun <- link$linkfun
  linkinv <- link$linkinv

#________________AR___________________________________
  if (isphi == T && istheta == F) {

    ar <- seq_along(phi)
    p <- max(ar)
    m <- 2*p

    ynew <- rep(alpha,(n+m))
    mu <- linkinv(ynew)

    eta <- y <- rep(NA, n+m)

    for(i in (m+1):(n+m)){
      eta[i]  <- alpha + (phi %*% ynew[i-ar])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- MxARMA::rmax(mu[i])
      ynew[i] <- linkfun(y[i])
    }
    return(y[(m+1):(n+m)])
  }

#_________________MA___________________________________________
  if (isphi == F && istheta == T){

    ma <- seq_along(theta)
    q <- max(ma)
    m <- 2*q

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    eta <- y <- error <- rep(0,n+m)

    for(i in (m+1):(n+m)){
      eta[i]  <- alpha + (theta %*% error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- MxARMA::rmax(mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]

    }
    return(y[(m+1):(n+m)])

  }

#________________ARMA________________________________________
  if(isphi == T && istheta == T){

    ar <- seq_along(phi)
    ma <- seq_along(theta)
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)

    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)

    error <- rep(0,n+m)
    eta <- y <- rep(NA, n+m)


    for(i in (m+1):(n+m)){
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- MxARMA::rmax(mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
    }

    return(y[(m+1):(n+m)])
  }
}
