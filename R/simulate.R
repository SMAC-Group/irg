#' Simulate Bivariate Irregularly Sampled Signals
#' @description This function allows to simulate a bivariate time series sampled at irregular intervals and following a bivariate first-order autoregressive process as
#' described in Heerah et al. (2020).
#' @param theta A \code{double} value of the vector representing the parameter vector for the model (the vector is length 8). The first four parameters are those for the first signal (root) and
#' represent respectively: \eqn{\phi}{phi} (the range parameter), \eqn{\psi}{psi} (the intensity parameter of the causal impact of the second signal on the first),
#' \eqn{\gamma}{gamma} (the parameter representing time of maximal impact of the second signal on the first), \eqn{\sigma^2}{sigma2} (the variance parameter of the first signal).
#' The following four parameters of the vector \code{theta} are the same but for the second signal.
#' @param times A \code{double} value of the vector containing the irregular time points (in the same unit) at which measurements are to be simulated.
#' @param burn_length An \code{integer} value representing how many measurements to simulate for the burn-in phase before taking the following values as the actual simulated time series
#' (defaults to \code{3*length(times)}).
#' @param burn_step A \code{double} value representing the fixed regular time points to use for the burn-in phase of the simulation (defaults to \code{median(diff(times))}).
#' @return A \code{list} containing the following objects:
#' \describe{
#'  \item{root}{The simulated values for the first signal (root)}
#'  \item{shoot}{The simulated values for the second signal (shoot)}
#' }
#' @author Roberto Molinari and Stephane Guerrier
#' @export
#' @examples
#' theta <- c(1, 0.99, 10, 0.01, 1, 0, 0, 0.1)
#' times <- c(0, 5, 10, 15, 20, 30, 45, 60, 90, 120)
#' sim_proc(theta, times)
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
