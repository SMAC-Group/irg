% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{sim_proc}
\alias{sim_proc}
\title{Simulate Bivariate Irregularly Sampled Signals}
\usage{
sim_proc(
  theta,
  times,
  burn_length = 3 * length(times),
  burn_step = median(diff(times))
)
}
\arguments{
\item{theta}{A \code{double} value of the vector representing the parameter vector for the model (the vector is length 8). The first four parameters are those for the first signal (root) and
represent respectively: \eqn{\phi}{phi} (the range parameter), \eqn{\psi}{psi} (the intensity parameter of the causal impact of the second signal on the first),
\eqn{\gamma}{gamma} (the parameter representing time of maximal impact of the second signal on the first), \eqn{\sigma^2}{sigma2} (the variance parameter of the first signal).
The following four parameters of the vector \code{theta} are the same but for the second signal.}

\item{times}{A \code{double} value of the vector containing the irregular time points (in the same unit) at which measurements are to be simulated.}

\item{burn_length}{An \code{integer} value representing how many measurements to simulate for the burn-in phase before taking the following values as the actual simulated time series
(defaults to \code{3*length(times)}).}

\item{burn_step}{A \code{double} value representing the fixed regular time points to use for the burn-in phase of the simulation (defaults to \code{median(diff(times))}).}
}
\value{
A \code{list} containing the following objects:
\describe{
 \item{root}{The simulated values for the first signal (root)}
 \item{shoot}{The simulated values for the second signal (shoot)}
}
}
\description{
This function allows to simulate a bivariate time series sampled at irregular intervals and following a bivariate first-order autoregressive process as
described in Heerah et al. (2020).
}
\examples{
theta <- c(1, 0.99, 10, 0.01, 1, 0, 0, 0.1)
times <- c(0, 5, 10, 15, 20, 30, 45, 60, 90, 120)
sim_proc(theta, times)
}
\author{
Roberto Molinari and Stephane Guerrier
}
