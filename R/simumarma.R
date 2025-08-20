#'  Simulation from the Maxwell ARMA (MxARMA) Model
#'
#'  Generates a univariate time series from an ARMA-type model whose innovations follow
#' the Maxwell distribution. The user can simulate pure autoregressive (AR),
#' pure moving average (MA), or mixed ARMA dynamics depending on which parameters
#' are provided.
#'
#' @param n Integer. Number of observations to simulate.
#' @param alpha  Numeric. Intercept of the linear predictor.
#' @param phi Numeric vector. Autoregressive coefficients. If \code{NULL} (default),
#'   no AR structure is included.
#' @param theta Numeric vector. Moving average coefficients. If \code{NULL} (default),
#'   no MA structure is included.
#'
#' @details
#' The MxARMA(p,q) model assumes that the conditional mean of the process is linked
#' to a linear predictor via a logarithmic link function. The innovations are assumed
#' to follow a Maxwell distribution, generated internally with \code{rmax()}.
#'
#' When only \code{phi} is supplied, an MxAR(p) model is simulated. When only
#' \code{theta} is supplied, an MxMA(q) model is generated. When both are supplied,
#' an MxARMA(p,q) model is simulated.
#'
#' @return A numeric vector of length \code{n} containing the simulated time series.
#'
#' @author Arthur Hintz
#'
#' @seealso \code{\link{rmax}}
#'
#' @examples
#' # AR(2) Example
#' set.seed(123)
#' ts_ar <- mxarma.sim(n = 500, alpha = 0.6, phi = c(0.7, -0.2))
#'
#' # MA(1) Example
#' ts_ma <- mxarma.sim(n = 500, alpha = 0.6, theta = 0.4)
#'
#' # ARMA(2,1) Example
#' ts_arma <- mxarma.sim(n = 500, alpha = 0.6, phi = c(0.5, 0.2), theta = 0.3)
#'
#' @importFrom stats make.link
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
      y[i]    <- MxARMA::rmax(mu = mu[i])
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
      y[i]    <- MxARMA::rmax(mu = mu[i])
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
      y[i]    <- MxARMA::rmax(mu = mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
    }

    return(y[(m+1):(n+m)])
  }
}
