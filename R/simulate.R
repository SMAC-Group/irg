#' Simulate Bivariate Irregularly Sampled Signals
#' @description TO DO
#' @param theta A \code{double} value of the vector containing the detrended and standardized measurements of the first signal (the root expressions).
#' @param times A \code{double} value of the vector containing the detrended and standardized measurements of the second signal (the shoot expressions).
#' @param burn_lenght A \code{double} value of the vector collecting the times (measured in the same unit) at which the measurements were taken
#' for both signals (\code{root} and \code{shoot}).
#' @param burn_step A \code{double} value for the starting values of the parameter vector \eqn{\theta}{theta} of the full causal model (must be of length 8).
#' @note XXX
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{stat}{The value of the Granger-Causal test statistic for the specified alternative hypothesis \code{HA}}
#'  \item{parameters}{A \code{double} vector containing the values of the parameters relevant for the specified alternative hypothesis \code{HA}
#'  (i.e. intensity of the causal impact (\eqn{\psi}{psi}) and the time of maximal impact (\eqn{\gamma}{gamma}))}
#' }
#' @details
#' XXX
#' @author Roberto Molinari and Stephane Guerrier
#' @export
sim_proc <- function(theta, times, burn_length = 3*length(times), burn_step = median(diff(times))) {

  phi.root <- theta[1]
  psi.root <- theta[2]
  delta.root <- theta[3]
  sigma.root <- theta[4]
  phi.shoot <- theta[5]
  psi.shoot <- theta[6]
  delta.shoot <- theta[7]
  sigma.shoot <- theta[8]

  shoot.star <- rnorm(1)
  root.star <- rnorm(1)

  del.t <- diff(times)
  del.t2 <- c(rep(burn_step, burn_length), del.t)

  for(i in 2:length(del.t2)) {

    shoot.star <- c(shoot.star, phi_func(del.t2[i], phi.shoot)*shoot.star[i-1] + psi_func(del.t2[i], delta.shoot, psi.shoot)*root.star[i-1] + rnorm(1, 0, sigma.shoot*(1-exp(-2*del.t2[i]/phi.shoot))))
    root.star <- c(root.star, phi_func(del.t2[i], phi.root)*root.star[i-1] + psi_func(del.t2[i], delta.root, psi.root)*shoot.star[i-1] + rnorm(1, 0, sigma.root*(1-exp(-2*del.t2[i]/phi.root))))

  }

  start_index <- length(del.t2) - length(del.t)
  root <- root.star[start_index:length(del.t2)]
  shoot <- shoot.star[start_index:length(del.t2)]

  out <- list(root = root, shoot = shoot)

  return(out)

}
