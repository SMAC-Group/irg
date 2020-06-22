# The f() function
phi_func <- function(del.t, phi) {

  return( exp(-del.t/phi) )

}

# The g() function
psi_func <- function(del.t, delta, psi) {

  eta <- max(del.t)

  return( psi*exp(-(del.t - delta)^2/eta) )

}


# The g() function
# psi_func <- function(del.t, delta, psi, eta) {
#
#   return( psi*exp(-(del.t - delta)^2/eta) )
#
# }

# The mean function for the likelihood
mu_func <- function(del.t, phi, delta, psi, x, y) {

  return( phi_func(del.t, phi)*x + psi_func(del.t, delta, psi)*y )

}

# The variance (sigma) function for the likelihood
sig_func <- function(del.t, phi, sigma) {

  return( sigma*(1 - exp(-2*del.t/phi)) )

}
