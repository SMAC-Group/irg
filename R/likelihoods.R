# Likelihood under null hypothesis (no impact)
lik_nodir <- function(theta, root, shoot, del.t) {

  phi.root <- exp(theta[1])
  psi.root <- 0
  delta.root <- 0
  sigma.root <- exp(theta[2])
  phi.shoot <- exp(theta[3])
  psi.shoot <- 0
  delta.shoot <- 0
  sigma.shoot <- exp(theta[4])

  iter1 <- 0

  for (i in 2:length(root)) {

    mu_t <- mu_func(del.t[i-1], phi.root, delta.root, psi.root, root[i-1], shoot[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.root, sigma.root)
    fi <- log(sig2_t) + (root[i] - mu_t)^2/sig2_t
    iter1 <- iter1 + fi

  }

  iter2 <- 0

  for (i in 2:length(root)) {

    mu_t <- mu_func(del.t[i-1], phi.shoot, delta.shoot, psi.shoot, shoot[i-1], root[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.shoot, sigma.shoot)
    fi <-  log(sig2_t) + (shoot[i] - mu_t)^2/sig2_t
    iter2 <- iter2 + fi

  }

  iter <- iter1 + iter2

  return(iter)

}


# Likelihood for both directions
lik_twodir <- function(theta, root, shoot, del.t){

  phi.root <- exp(theta[1])
  psi.root <- logit2_inv(theta[2])/2
  delta.root <- exp(theta[3])
  sigma.root <- exp(theta[4])
  phi.shoot <- exp(theta[5])
  psi.shoot <- logit2_inv(theta[6])/2
  delta.shoot <- exp(theta[7])
  sigma.shoot <- exp(theta[8])

  # phi.root <- exp(theta[1])
  # psi.root <- logit2_inv(theta[2])/2
  # delta.root <- 55*logit_inv(theta[3]) + 5
  # sigma.root <- exp(theta[4])
  # phi.shoot <- exp(theta[5])
  # psi.shoot <- logit2_inv(theta[6])/2
  # delta.shoot <- 55*logit_inv(theta[7]) + 5
  # sigma.shoot <- exp(theta[8])

  iter1 <- 0

  for (i in 2:length(root)) {

    mu_t <- mu_func(del.t[i-1], phi.root, delta.root, psi.root, root[i-1], shoot[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.root, sigma.root)
    fi <- log(sig2_t) + (root[i] - mu_t)^2/sig2_t
    iter1 <- iter1 + fi

  }

  iter2 <- 0

  for (i in 2:length(root)) {

    mu_t <- mu_func(del.t[i-1], phi.shoot, delta.shoot, psi.shoot, shoot[i-1], root[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.shoot, sigma.shoot)
    fi <- log(sig2_t) + (shoot[i] - mu_t)^2/sig2_t
    iter2 <- iter2 + fi

  }

  iter <- iter1 + iter2

  return(iter)

}


# Likelihood for shoot to root
lik_shoot_root <- function(theta, root, shoot, del.t){

  phi.root <- exp(theta[1])
  psi.root <- logit2_inv(theta[2])/2
  delta.root <- exp(theta[3])
  sigma.root <- exp(theta[4])
  phi.shoot <- exp(theta[5])
  psi.shoot <- 0
  delta.shoot <- 0
  sigma.shoot <- exp(theta[6])

  # phi.root <- exp(theta[1])
  # psi.root <- logit2_inv(theta[2])/2
  # delta.root <- (55)*logit_inv(theta[3]) + 5
  # sigma.root <- exp(theta[4])
  # phi.shoot <- exp(theta[5])
  # psi.shoot <-  0
  # delta.shoot <- 10
  # sigma.shoot <- exp(theta[6])

  iter1 <- 0

  for (i in 2:length(root)) {

    mu_t <- mu_func(del.t[i-1], phi.root, delta.root, psi.root, root[i-1], shoot[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.root, sigma.root)
    fi <- log(sig2_t) + (root[i] - mu_t)^2/sig2_t
    iter1 <- iter1 + fi

  }

  iter2 <- 0

  for (i in 2:length(root)) {

    mu_t <- mu_func(del.t[i-1], phi.shoot, delta.shoot, psi.shoot, shoot[i-1], root[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.shoot, sigma.shoot)
    fi <- log(sig2_t) + (shoot[i] - mu_t)^2/sig2_t
    iter2 <- iter2 + fi

  }

  iter <- iter1 + iter2

  return(iter)

}


# Likelihood for root to shoot
lik_root_shoot <- function(theta, root, shoot, del.t) {

  phi.root <- exp(theta[1])
  psi.root <- 0
  delta.root <- 0
  sigma.root <- exp(theta[2])
  phi.shoot <- exp(theta[3])
  psi.shoot <- logit2_inv(theta[4])/2
  delta.shoot <- exp(theta[5])
  sigma.shoot <- exp(theta[6])

  # phi.root <- exp(theta[1])
  # psi.root <- 0
  # delta.root <- 10
  # sigma.root <- exp(theta[2])
  # phi.shoot <- exp(theta[3])
  # psi.shoot <- logit2_inv(theta[4])/2
  # delta.shoot <- 55*logit_inv(theta[5]) + 5
  # sigma.shoot <- exp(theta[6])

  iter1 <- 0

  for (i in 2:length(root)){

    mu_t <- mu_func(del.t[i-1], phi.root, delta.root, psi.root, root[i-1], shoot[i-1])
    sig2_t <- sig_func(del.t[i-1], phi.root, sigma.root)
    fi <- log(sig2_t) + (root[i] - mu_t)^2/sig2_t
    iter1 <- iter1 + fi

  }

  iter2 <- 0

  for (i in 2:length(root)){

    mu_t <- mu_func(del.t[i-1], phi.shoot, delta.shoot, psi.shoot, shoot[i-1], root[i-1])
    sig2_t <- sig_func(del.t[i-1],phi.shoot,sigma.shoot)
    fi <- log(sig2_t) + (shoot[i] - mu_t)^2/sig2_t
    iter2 <- iter2 + fi

  }

  iter <- iter1 + iter2

  return(iter)

}
