#' Probability Density Function of the Parametrized Maxwell Distribution
#'
#' Computes the probability density function (PDF) of a random variable \eqn{Y}
#' that follows a Maxwell distribution with a scale parameter \eqn{\mu > 0}.
#' The density function is given by:
#'
#' \deqn{f(y;\mu) = \frac{32y^2}{\pi^2\mu^3} \exp\left(-\frac{4y^2}{\pi\mu^2}\right)}
#'
#' This distribution is defined for \eqn{y \ge 0}.
#'
#' @param y A numeric value or vector at which the density is to be evaluated. Must be non-negative.
#' @param mu A positive numeric value representing the scale parameter \eqn{\mu}.
#'
#' @return A numeric vector with the values of the density function.
#'
#' @examples
#' dmax(y = 1, mu = 2)
#'
#' # Vectorized example
#' dmax(y = c(0.5, 1, 2), mu = 2)
#'
#' # Plotting the PDF for mu = 3
#' curve(dmax(x, 3), from = 0, to = 10,
#'       xlab = "y", ylab = "Density", n = 1000)
#'
#' @export
dmax <- function(y, mu) {
  if (!is.numeric(y)) stop("'y' must be numeric.")
  if (!is.numeric(mu) || length(mu) != 1 || mu <= 0)
    stop("'mu' must be a single positive number.")

  density <- ifelse(y < 0, 0,
                    (32 * y^2 / (pi^2 * mu^3)) * exp(-4 * y^2 / (pi * mu^2))
  )

  return(density)
}
