# Utility functions for hmm
# The "## ---- function_name" text lets us define code chunks to import and display in the supplement on GitHub easily
# without needing to copy-paste.

## ---- delta.n2w
# Function to transform natural parameters to working ones
delta.n2w <- function(m, delta) {
  tdelta <- log(delta[- 1] / delta[1])
  return(tdelta) 
}

## ---- delta.w2n
# Function to transform working parameters to natural ones
delta.w2n <- function(m, tdelta) {
  if (m == 1) return(1)
  
  # set first element to one and fill in the last m - 1 elements with working parameters and take exp
  foo <- c(1, exp(tdelta))
  
  # normalize
  delta <- foo / sum(foo)
  
  return(delta)
}

## ---- gamma.n2w
# Function to transform natural parameters to working ones
gamma.n2w <- function(m, gamma) {
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  return(tgamma)
}

## ---- gamma.w2n
# Function to transform working parameters to natural ones
gamma.w2n <- function(m, tgamma) {
  gamma <- diag(m)
  if (m == 1) return(gamma)
  gamma[!gamma] <- exp(tgamma)
  gamma <- gamma/apply(gamma, 1, sum)
  return(gamma)
}

## ---- quantile.colwise
# 2.5% and 97.5% quantiles
quantile.colwise <- function(data) {
  return(quantile(data, probs = c(0.05 / 2, 1 - 0.05 / 2)))
}

# Transform tpm with column wise idx to row and column indices
matrix.col.idx.to.rowcol <- function(idx, m) {
  # The indices are 1:m in the 1st column, then (m+1):(2*m) in the 2nd, etc...
  row <- (idx - 1) %% m + 1
  col <- (idx - 1) %/% m + 1
  return(c(row, col))
}

## ---- get.norm.emission.probs
# Calculate emission probabilities
get.norm.emission.probs <- function(data, mu, sigma) {
  n <- length(data)
  m <- length(mu)
  emission_probs <- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) {
    if (is.na(data[i])) {
      emission_probs[i, ] <- rep(1, m)
    } else {
      emission_probs[i, ] <- dnorm(data[i], mu, sigma)
    }
  }
  return(emission_probs)
}

## ---- stat.dist
# Compute the stationary distribution of a Markov chain
# with transition probability gamma
stat.dist <- function(gamma) {
  # The code from Zucchini can crash when dealing with computationally small numbers.
  # This is likely an approximation error.
  # For some reason, the result may be a complex number with imaginary part equal to 0.
  # We can just ignore the imaginary part because the eigenvectors of a real eigenvalue must be real.
  # In the cases where it happened, solve(t(diag(m) - gamma + 1), rep(1, m)) produced the same result without the imaginary part.
  first_eigen_row <- Re(eigen(t(gamma))$vectors[, 1])
  return(first_eigen_row / sum(first_eigen_row))
}

## ---- norm.HMM.pn2pw
# Transform Gaussian natural parameters to working parameters
norm.HMM.pn2pw <- function(m, mu, sigma, gamma, delta = NULL,
                           stationary = TRUE) {
  tsigma <- log(sigma)
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if (stationary) {
    # If tdelta is set to NULL and returned in the list,
    # it will cause issues when optimizing with TMB
    return(list(tmu = mu, tsigma = tsigma, tgamma = tgamma))
  } else {
    tdelta <- log(delta[- 1] / delta[1])
    # TMB requires a list
    return(list(tmu = mu, tsigma = tsigma, tgamma = tgamma, tdelta = tdelta))
  }
}

## ---- norm.HMM.pw2pn
# Transform Gaussian working parameters to natural parameters
norm.HMM.pw2pn <- function(m, parvect, stationary = TRUE) {
  mu_indices <- 1:m
  sigma_indices <- max(mu_indices) + 1:m
  gamma_indices <- max(sigma_indices) + 1:(m ^ 2)
  delta <- max(gamma_indices) + 1:m
  
  parvect <- unlist(parvect)
  mu <- parvect[mu_indices]
  sigma <- exp(parvect[sigma_indices])
  gamma <- diag(m)
  if (m == 1) return(list(tmu = mu, sigma = sigma, gamma = gamma, delta = 1))
  gamma[!gamma] <- exp(parvect[gamma_indices])
  gamma <- gamma / apply(gamma, 1, sum)
  if (stationary) {
    # The code from Zucchini can crash when dealing with computationally small numbers.
    delta <- stat.dist(gamma)
  } else {
    foo <- c(1, exp(parvect[(m * m + 1):(m * m + m - 1)]))
    delta <- foo / sum(foo)
  }
  return(list(tmu = mu, sigma = sigma, gamma = gamma, delta = delta))
}

## ---- norm.HMM.generate.sample
# Generate a random sample from a HMM
norm.HMM.generate.sample  <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  }
  x <- rnorm(ns, mean = mod$mu[state], sd = mod$sigma[state])
  return(list(data = x, state = state))
}

## ---- norm.HMM.label.order
# Relabel states by increasing Poisson means
norm.HMM.label.order <- function(m, mu, sigma, gamma, delta = NULL) {
  # Get the indexes of the sorted states
  # according to ascending lambda
  sorted_mu <- sort(mu, index.return = TRUE)$ix
  # Re-order the TPM according to the switched states
  # in the sorted lambda
  ordered_gamma <- matrix(0, nrow = m, ncol = m)
  for (col in 1:m) {
    new_col <- which(sorted_mu == col)
    for (row in 1:m) {
      new_row <- which(sorted_mu == row)
      ordered_gamma[row, col] <- gamma[new_row, new_col]
    }
  }
  ordered_sigma <- sigma[sorted_mu]
  # Re-order the stationary distribution if it was provided
  # Generate it otherwise
  if (is.null(delta)) {
    delta <- stat.dist(ordered_gamma)
  } else {
    delta <- delta[sorted_mu]
  }
  return(list(mu = sort(mu),
              sigma = ordered_sigma,
              gamma = ordered_gamma,
              delta = delta))
}
