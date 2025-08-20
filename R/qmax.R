#' Quantile Function of the Parametrized Maxwell Distribution
#'
#' Computes the quantile function of the Maxwell distribution with a scale parameter \eqn{\mu > 0}.
#' It returns the value \eqn{y} such that \eqn{P(Y \le y) = u}.
#'
#' @param u A numeric value or vector of probabilities in \eqn{[0, 1]}.
#' @param mu A positive numeric value representing the scale parameter \eqn{\mu}.
#'
#' @return A numeric vector of quantiles corresponding to the given probabilities.
#'
#' @examples
#' qmax(u = 0.5, mu = 2)
#' qmax(u = c(0.25, 0.5, 0.75), mu = 2)
#'
#' @export
qmax <- function(u, mu) {
  if (!is.numeric(u) || any(u < 0 | u > 1))
    stop("'u' must be numeric and between 0 and 1.")
  if (!is.numeric(mu) || length(mu) != 1 || mu <= 0)
    stop("'mu' must be a single positive number.")

  ext <- uniroot(function(x) MxARMA::pmax(x, mu) - u,
                 interval = c(0, mu * 20), tol = 1e-30,
                 extendInt = "yes")$root
  return(abs(ext))



}

