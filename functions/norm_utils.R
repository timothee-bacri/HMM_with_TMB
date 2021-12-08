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

## ---- TMB.estimate
# Estimation using TMB
TMB.norm.estimate <- function(TMB_data,
                              parameters,
                              MakeADFun_obj = NULL,
                              map = list(),
                              gradient = FALSE,
                              hessian = FALSE,
                              std_error = FALSE) {
  
  obj <- MakeADFun_obj
  if (is.null(MakeADFun_obj)) {
    obj <- MakeADFun(TMB_data, parameters, DLL = "norm_hmm", silent = TRUE, map = map)
  }
  
  # The function ifelse cannot return a NULL value, the function switch can
  # If gradient is FALSE, then gradient + 1 is 1 and the switch returns NULL
  # If gradient is TRUE, then gradient + 1 is 2 and the switch returns mod$gr
  gr <- switch(gradient + 1, NULL, obj$gr)
  he <- switch(hessian + 1, NULL, obj$he)
  
  m <- TMB_data$m
  
  # Optimizing while handling errors
  mod <- tryCatch({
    nlminb(start = obj$par, objective = obj$fn, gradient = gr, hessian = he)
  },
  error = function(e) {
    message("nlminb error:")
    message("m = ", m)
    message("gradient = ", gradient)
    message("hessian = ", hessian)
    message("The original error message is:\n", e)
    return()
  })
  convergence <- mod$convergence
  if (is.null(mod)) {
    return()
  }
  if (convergence != 0) {
    w <- paste0("nlminb didn't succesfully converge:\n",
                mod$message,
                "\nm = ", m,
                "\ngradient = ", gradient,
                "\nhessian = ", hessian, "\n")
    warning(w)
  }
  
  mllk <- mod$objective
  np <- length(unlist(parameters))
  AIC <- 2 * (mllk + np)
  n <- sum(!is.na(TMB_data$x))
  BIC <- 2 * mllk + np * log(n)
  # Return standard errors
  if (std_error) {
    adrep <- summary(sdreport(obj), "report")
    
    rows <- rownames(adrep) == "mu"
    mu <- adrep[rows, "Estimate"]
    mu_std_error <- adrep[rows, "Std. Error"]
    
    rows <- rownames(adrep) == "sigma"
    sigma <- adrep[rows, "Estimate"]
    sigma_std_error <- adrep[rows, "Std. Error"]
    
    rows <- rownames(adrep) == "gamma"
    gamma <- adrep[rows, "Estimate"]
    gamma <- matrix(gamma, ncol = m)
    gamma_std_error <- adrep[rows, "Std. Error"]
    gamma_std_error <- matrix(gamma_std_error, ncol = m)
    
    rows <- rownames(adrep) == "delta"
    delta <- adrep[rows, "Estimate"]
    delta_std_error <- adrep[rows, "Std. Error"]
    
    return(list(m = m, mu = mu, sigma = sigma, gamma = gamma, delta = delta,
                mu_std_error = mu_std_error,
                sigma_std_error = sigma_std_error,
                gamma_std_error = gamma_std_error,
                delta_std_error = delta_std_error,
                convergence = convergence, mllk = mllk,
                AIC = AIC, BIC = BIC,
                mod = mod, obj = obj))
  }
  
  mu <- as.numeric(mod$par[names(mod$par) == "mu"])
  tsigma <- as.numeric(mod$par[names(mod$par) == "tsigma"])
  tgamma <- as.numeric(mod$par[names(mod$par) == "tgamma"])
  
  estim <- norm.HMM.pw2pn(m, c(mu, tsigma, tgamma))
  
  return(list(m = m, mu = estim$mu, sigma = estim$sigma, gamma = estim$gamma,
              delta = estim$delta, convergence = convergence, mllk = mllk,
              AIC = AIC, BIC = BIC, mod = mod, obj = obj))
}

## ---- HMM.decode
# Computes log-forward, log-backward and conditional probabilities
# and decoding based on (optimized) MakeADFun object,
HMM.decode <- function(obj) {
  
  # Setup
  # Retrieve the objects at ML value
  adrep <- obj$report(obj$env$last.par.best)
  delta <- adrep$delta
  gamma <- adrep$gamma
  emission_probs <- adrep$emission_probs
  n <- adrep$n
  m <- length(delta)
  mllk <- adrep$mllk
  
  # Compute log-forward probabilities (scaling used)
  lalpha <- matrix(NA, m, n)
  foo <- delta * emission_probs[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- log(foo) + lscale
  for (i in 2:n) {
    foo <- foo %*% gamma * emission_probs[i, ]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  # Compute log-backwards probabilities (scaling used)
  lbeta <- matrix(NA, m, n)
  lbeta[, n] <- rep(0, m)
  foo <- rep (1 / m, m)
  lscale <- log(m)
  for (i in (n - 1):1) {
    foo <- gamma %*% (emission_probs[i + 1, ] * foo)
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  # Compute conditional state probabilities, smoothing probabilities
  stateprobs <- matrix(NA, ncol = n, nrow = m)
  llk <- - mllk
  for(i in 1:n) {
    stateprobs[, i] <- exp(lalpha[, i] + lbeta[, i] - llk)
  }
  
  # Most probable states
  ldecode <- rep(NA, n)
  for (i in 1:n) {
    ldecode[i] <- which.max(stateprobs[, i])
  }
  
  # Output as list
  list(lalpha = lalpha, lbeta = lbeta,
       stateprobs = stateprobs, ldecode = ldecode)
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

## ---- norm.HMM.mllk
# Calculate the negative log-likelihood, based on the book
norm.HMM.mllk <- function(parvect, x_alias, m_alias, stationary = TRUE) {
  # The variable names m and x are already used as parameters for the hessian
  # m_alias and x_alias are only replacement names
  m <- m_alias
  x <- x_alias
  n <- length(x)
  pn <- norm.HMM.pw2pn(m, parvect)
  emission_probs <- get.norm.emission.probs(x, pn$mu, pn$sigma)
  
  if (m == 1) return(- sum(log(emission_probs[, 1])))
  
  foo <- pn$delta * emission_probs[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  for (i in 2:n) {
    if (!is.na(x[i])) {
      P <- emission_probs[i, ]
    } else {
      P <- rep(1, m)
    }
    
    foo <- foo %*% pn$gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  mllk <- - lscale
  return(mllk)
}

## ---- DM.estimate
# Compute the ML estimates without using TMB
DM.norm.estimate = function(x, m, mu0, sigma0, gamma0, delta0 = NULL,
                            stationary = TRUE) {
  parvect0 <- norm.HMM.pn2pw(m = m, mu = mu0, sigma = sigma0, gamma = gamma0,
                             delta = delta0, stationary = stationary)
  # nlminb needs a vector, not a list
  parvect0 <- unlist(parvect0)
  
  mod <- nlminb(start = parvect0, objective = norm.HMM.mllk, x_alias = x, m_alias = m, stationary = stationary)
  pw <- mod$par
  pn <- norm.HMM.pw2pn(m, as.numeric(pw))
  mllk <- mod$objective
  np <- unlist(parvect0)
  np <- length(np)
  AIC <- 2 * (mllk + np)
  n <- sum(!is.na(x))
  BIC <- 2 * mllk + np * log(n)
  return(list(m = m, mu = pn$mu, sigma = pn$sigma, gamma = pn$gamma, delta = pn$delta,
              convergence = mod$convergence, mllk = mllk, AIC = AIC, BIC = BIC,
              mod = mod))
  
}
## ---- nlmfn
# Function to use TMB's gradient and/or hessian in nlm
nlmfn <- function(par, obj, gr = TRUE, he = TRUE) {
  res <- as.numeric(obj$fn(par))
  if(gr) {attr(res, "gradient") <- obj$gr(par)}
  if(he) {attr(res, "hessian") <- obj$he(par)}
  res
}

