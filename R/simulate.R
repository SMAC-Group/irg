# Simulate process from model
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
