#' Cumulative Distribution Function of the Parametrized Maxwell Distribution
#'
#' Computes the cumulative distribution function (CDF) of a random variable \eqn{Y}
#' that follows a Maxwell distribution with a scale parameter \eqn{\mu > 0}.
#' The CDF is given by:
#'
#' \deqn{F(y;\mu) = \frac{2}{\sqrt{\pi}} \cdot \gamma\left(\frac{3}{2}, \frac{4y^2}{\pi \mu^2} \right)}
#' where \eqn{\gamma(a, x)} is the lower incomplete gamma function (regularized via \code{pgamma}).
#'
#' This distribution is defined for \eqn{y \ge 0}.
#'
#' @param y A numeric value or vector at which the CDF is to be evaluated. Must be non-negative.
#' @param mu A positive numeric value representing the scale parameter \eqn{\mu}.
#'
#' @return A numeric vector of cumulative probabilities.
#'
#' @examples
#' pmax(y = 1, mu = 2)
#'
#' # Vectorized input
#' pmax(y = c(0.1, 0.5, 1), mu = 2)
#'
#' # Plotting the CDF for mu = 3
#' curve(pmax(x, 3), from = 0, to = 10,
#'       xlab = "y", ylab = "CDF", n = 1000)
#'
#' @export
pmax <- function(y, mu) {
  if (!is.numeric(y)) stop("'y' must be numeric.")
  if (!is.numeric(mu) || length(mu) != 1 || mu <= 0)
    stop("'mu' must be a single positive number.")

  cdf <- ifelse(y < 0, 0,
                (2/sqrt(pi)) * gamma(3/2) * pgamma((4 * y^2)/(pi * mu^2),
                                                   3/2, lower.tail = TRUE))

  return(cdf)
}

