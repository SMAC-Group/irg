# Functions to estimate starting values
ar1_to_gm <- function(theta, freq) {

  beta <- -log(theta[1])*freq
  sigma2_gm <- theta[2]/(1 - theta[1]^2)

  return(c(beta, sigma2_gm))

}

gm_to_ar1 <- function(theta, freq) {

  phi <- exp(-theta[1]/freq)
  sigma2 <- theta[2]*(1 - exp(-2*theta[1]/freq))

  return(c(phi, sigma2))

}

logit <- function(x) {

  log(x/(1 - x))

}

logit_inv <- function(x) {

  exp(x)/(exp(x) + 1)

}


logit2 <- function(x) {

  log((x + 2)/(2 - x))

}

logit2_inv <- function(x) {

  4/(1 + exp(-x)) - 2

}

start_values <- function(root, shoot, del.t) {

  fit.root <- lm(root[2:length(root)] ~ root[1:(length(root)-1)] + shoot[1:(length(shoot)-1)] - 1)
  fit.shoot <- lm(shoot[2:length(shoot)] ~ shoot[1:(length(shoot)-1)] + root[1:(length(root)-1)] - 1)
  phi.gm.root <- abs(ar1_to_gm(c(abs(fit.root$coefficients[1]), (summary(fit.root)$sigma)**2), median(del.t)))
  psi.root <- fit.root$coefficients[2]
  phi.gm.shoot <- abs(ar1_to_gm(c(abs(fit.shoot$coefficients[1]), (summary(fit.shoot)$sigma)**2), median(del.t)))
  psi.shoot <- fit.shoot$coefficients[2]

  theta.start <- as.numeric(c(log(phi.gm.root[1]), logit2(ifelse(abs(psi.root) >= 1, sign(psi.root)*0.99, psi.root)*2), log(median(del.t)), log(phi.gm.root[2]),
                              log(phi.gm.shoot[1]), logit2(ifelse(abs(psi.shoot) >= 1, sign(psi.shoot)*0.99, psi.shoot)*2), log(median(del.t)), log(phi.gm.shoot[2])))

  return(theta.start)

}
