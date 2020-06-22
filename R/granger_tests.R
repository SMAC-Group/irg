#' Granger-Causal Test Statistic
#' @description This function is used within the \code{granger_test} function to compute the Granger-Causal test statistic on the signals of interest
#' as well as on the simulated ones based on the \code{sim_proc} function (see Heerah et al., 2020).
#' @param root A \code{double} value of the vector containing the detrended and standardized measurements of the first signal (the root expressions).
#' @param shoot A \code{double} value of the vector containing the detrended and standardized measurements of the second signal (the shoot expressions).
#' @param times A \code{double} value of the vector collecting the times (measured in the same unit) at which the measurements were taken
#' for both signals (\code{root} and \code{shoot}).
#' @param theta A \code{double} value for the starting values of the parameter vector \eqn{\theta}{theta} of the full causal model (must be of length 8).
#' @param alternative A \code{character} value defining the alternative hypothesis for the Granger-Causal test (with the null hypothesis always being "the signals do not Granger-cause each other").
#' The default value is \code{"twodir"} which tests the alternative that "both signals Granger-cause each other". The other options are \code{"rtos"} which tests the alternative
#' "the first signal (root) Granger-causes the second signal (shoot)" and \code{"stor"} which tests the alternative "the second signal (shoot) Granger-causes the first signal (root)".
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
causal_stat <- function(root, shoot, times, theta, alternative = "twodir") {

  del.t <- diff(times)

  h0 <- optim(theta[c(1, 4, 5, 8)], lik_nodir, root = root, shoot = shoot, del.t = del.t)$value

  if(alternative == "twodir") {

    hA <- optim(theta, lik_twodir, root = root, shoot = shoot, del.t = del.t)

    granger_param <- c(logit2_inv(hA$par[2])/2, exp(hA$par[3]), logit2_inv(hA$par[6])/2, exp(hA$par[7]))
    names(granger_param) <- c("Shoot to Root Intensity", "Shoot to Root Impact Time", "Root to Shoot Intensity", "Root to Shoot Impact Time")

  } else {

    if(alternative == "rtos") {

      hA <- optim(theta[c(1, 4, 5:8)], lik_root_shoot, root = root, shoot = shoot, del.t = del.t)

      granger_param <- c(logit2_inv(hA$par[4])/2, exp(hA$par[5]))
      names(granger_param) <- c("Root to Shoot Intensity", "Root to Shoot Impact Time")

    } else {

      if(alternative == "stor") {

        hA <- optim(theta[c(1:4, 5, 8)], lik_shoot_root, root = root, shoot = shoot, del.t = del.t)

        granger_param <- c(logit2_inv(hA$par[2])/2, exp(hA$par[3]))
        names(granger_param) <- c("Shoot to Root Intensity", "Shoot to Root Impact Time")

      } else {

        stop("Alternative hypothesis is not well specified")

      }

    }

  }

  granger_stat <- h0 - hA$value

  return(list("stat" = granger_stat, "parameters" = granger_param))

}


#' @title Granger-Causal Test for Irregularly Sampled Signals
#' @description This function takes two signals collected at the same irregularly spaced intervals and performs a Granger-Causal test for the required alternative hypothesis
#' (i.e. one signal Granger-causes the other or both signals Granger-cause each other). See Heerah et al., 2020 for more details on the procedure.
#' @param root A \code{double} value of the vector containing the detrended and standardized measurements of the first signal (the root expressions).
#' @param shoot A \code{double} value of the vector containing the detrended and standardized measurements of the second signal (the shoot expressions).
#' @param times A \code{double} value of the vector collecting the times (measured in the same unit) at which the measurements were taken
#' for both signals (\code{root} and \code{shoot}).
#' @param theta A \code{double} value of the vector providing the starting values for the full model parameter vector (for this reason the length of the vector must be exactly 8).
#' If not provided (default is \code{theta = NULL}), then starting values will be found using an approximation based on the assumption of regularly spaced measurements.
#' @param H A positive \code{integer} representing the number of parametric bootstrap replicates used to estimate the distribution of the test statistic under the null hypothesis.
#' @param alternative A \code{character} value defining the alternative hypothesis for the Granger-Causal test (with the null hypothesis always being "the signals do not Granger-cause each other").
#' The default value is \code{"twodir"} which tests the alternative that "both signals Granger-cause each other". The other options are \code{"rtos"} which tests the alternative
#' "the first signal (root) Granger-causes the second signal (shoot)" and \code{"stor"} which tests the alternative "the second signal (shoot) Granger-causes the first signal (root)".
#' @param seed A positive \code{integer} defining the seed value for the simulated signals throughout the parametric bootstrap procedure (default is \code{seed = 123}).
#' @note XXX
#' @return A \code{list} with the following objects:
#' \describe{
#'  \item{alternative}{The code of the alternative hypothesis that is tested.}
#'  \item{pvalue}{The p-value computed using the bootstrap distribution.}
#'  \item{parameters}{The parameters of "intensity" and "time of impact" for the alternative hypothesis that is tested.}
#' }
#' @details
#' The procedure is based on the work in Heerah et al., 2020 where the null hypothesis states that "neither of the signals Granger-causes the other".
#' According to the alternative hypothesis specified by the user, the function computes the corresponding p-value using a parametric bootstrap procedure.
#' If the null hypothesis is rejected in favour of the alternative, then the user can obtain information on (i) how strong the impact of a signal on the other is and
#' (ii) its direction (positive or negative) also referred to as "intensity" (represented by \eqn{\psi}{psi}) as well as the time point at which this impact is maximal
#' also referred to as "impact time" (represented by \eqn{\gamma}{gamma}).
#' @author Roberto Molinari and Stephane Guerrier
#' @export
#' @examples
#' \dontrun{
#' data(signals)
#' times <- c(0, 5, 10, 15, 20, 30, 45, 60, 90, 120)
#' granger_test(root = signals$root, shoot = signals$shoot, times = times, alternative = "twodir")
#' }
granger_test <- function(root, shoot, times, theta = NULL, alternative = "twodir", H = 100, seed = 123) {

  del.t <- diff(times)

  if(is.numeric(theta)) {

    if(length(theta) != 8) stop("Initial parameter vector must contain exactly 8 numeric values")

    theta_start <- c(log(theta[1]), logit2(theta[2]*2), log(theta[3]), log(theta[4]), log(theta[5]), logit2(theta[6]*2), log(theta[7]), log(theta[8]))

  } else {

    theta_start <- start_values(root = root, shoot = shoot, del.t = del.t)

  }

  # Get test statistic
  test_stat <- causal_stat(root, shoot, times, theta_start, alternative = alternative)

  # Estimate model under H0
  theta_hat <- optim(theta_start[c(1, 4, 5, 8)], lik_nodir, root = root, shoot = shoot, del.t = del.t)$par
  theta_h0 <- c(exp(theta_hat[1]), 0, 0, exp(theta_hat[2]), exp(theta_hat[3]), 0, 0, exp(theta_hat[4]))

  # Initialisation for parametric bootstrap
  boot_dist <- rep(NA, H)

  # Start bootstrap
  pb <- txtProgressBar(min = 0, max = H, style = 3)

  for (i in 1:H) {

    set.seed(i + seed)
    sim_boot <- sim_proc(theta_h0, times)
    boot_dist[i] <- causal_stat(sim_boot$root, sim_boot$shoot, times, theta_h0, alternative = alternative)$stat

    # update progress bar
    setTxtProgressBar(pb, i)

  }

  close(pb)

  # Refresh seed for system
  set.seed(Sys.time())

  # Compute p-value
  p_value <- (sum(test_stat$stat < boot_dist) + 1)/(H + 1)

  # Output results
  output <- list("alternative" = alternative, "pvalue" = p_value, "parameters" = test_stat$parameters)

  return(output)

}
