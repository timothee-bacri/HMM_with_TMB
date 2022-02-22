# Utility functions for hmm
# The "## ---- function_name" text lets us define code chunks to import and display in the supplement on GitHub easily
# without needing to copy-paste.

## ---- delta.n2w
# Function to transform natural parameters to working ones
delta.n2w <- function(m, delta){
  tdelta <- log(delta[- 1] / delta[1])
  return(tdelta) 
}

## ---- delta.w2n
# Function to transform working parameters to natural ones
delta.w2n <- function(m, tdelta){
  if (m == 1) return(1)
  
  # set first element to one and fill in the last m - 1 elements with working parameters and take exp
  foo <- c(1, exp(tdelta))
  
  # normalize
  delta <- foo / sum(foo)
  
  return(delta)
}

## ---- gamma.n2w
# Function to transform natural parameters to working ones
gamma.n2w <- function(m, gamma){
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  return(tgamma)
}

## ---- gamma.w2n
# Function to transform working parameters to natural ones
gamma.w2n <- function(m, tgamma){
  gamma <- diag(m)
  if (m == 1) return(gamma)
  gamma[!gamma] <- exp(tgamma)
  gamma <- gamma / apply(gamma, 1, sum)
  return(gamma)
}

## ---- TMB.estimate
# Estimation using TMB, wrapper of stats::nlminb capable of moderate error handling
TMB.estimate <- function(TMB_data,
                         parameters,
                         map = list(),
                         gradient = FALSE,
                         hessian = FALSE,
                         std_error = FALSE) {
  
  obj <- MakeADFun(TMB_data, parameters, DLL = "poi_hmm", silent = TRUE, map = map)
  
  gr <- if (gradient) obj$gr else NULL
  he <- if (hessian) obj$he else NULL
  
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
    
    rows <- rownames(adrep) == "lambda"
    lambda <- adrep[rows, "Estimate"]
    lambda_std_error <- adrep[rows, "Std. Error"]
    
    rows <- rownames(adrep) == "gamma"
    gamma <- adrep[rows, "Estimate"]
    gamma <- matrix(gamma, ncol = m)
    gamma_std_error <- adrep[rows, "Std. Error"]
    gamma_std_error <- matrix(gamma_std_error, ncol = m)
    
    rows <- rownames(adrep) == "delta"
    delta <- adrep[rows, "Estimate"]
    delta_std_error <- adrep[rows, "Std. Error"]
    
    return(list(m = m, lambda = lambda, gamma = gamma, delta = delta,
                lambda_std_error = lambda_std_error,
                gamma_std_error = gamma_std_error,
                delta_std_error = delta_std_error,
                convergence = convergence, mllk = mllk,
                AIC = AIC, BIC = BIC,
                mod = mod, obj = obj))
  }
  
  tlambda <- as.numeric(mod$par[names(mod$par) == "tlambda"])
  tgamma <- as.numeric(mod$par[names(mod$par) == "tgamma"])
  
  estim <- pois.HMM.pw2pn(m, c(tlambda, tgamma))
  
  return(list(m = m, lambda = estim$lambda, gamma = estim$gamma,
              delta = estim$delta, convergence = convergence, mllk = mllk,
              AIC = AIC, BIC = BIC, mod = mod, obj = obj))
}

## ---- HMM.decode
# Computes log-forward, log-backward and conditional probabilities
# and decoding based on an optimized TMB::MakeADFun object
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
# Computes the 2.5% and 97.5% quantiles
quantile.colwise <- function(data) {
  return(quantile(data, probs = c(0.05 / 2, 1 - 0.05 / 2)))
}

# Transform tpm with column wise idx to row and column indices
# When a matrix is passed where a vector is expected (for example in a column of a data.frame), the matrix is read column-wise.
matrix.col.idx.to.rowcol <- function(idx, m) {
  # The indices are 1:m in the 1st column, then (m+1):(2*m) in the 2nd, etc...
  row <- (idx - 1) %% m + 1
  col <- (idx - 1) %/% m + 1
  return(c(row, col))
}

## ---- get.emission.probs
# Calculate emission probabilities
get.emission.probs <- function(data, lambda) {
  n <- length(data)
  m <- length(lambda)
  emission_probs <- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) {
    if (is.na(data[i])) {
      emission_probs[i, ] <- rep(1, m)
    } else {
      emission_probs[i, ] <- dpois(data[i], lambda)
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
  # In the cases where it happened, solve(t(diag(m) - gamma + 1), rep(1, m)) produced the same result
  # without the imaginary part.
  first_eigen_row <- Re(eigen(t(gamma))$vectors[, 1])
  return(first_eigen_row / sum(first_eigen_row))
}

## ---- pois.HMM.pn2pw
# Transform Poisson natural parameters to working parameters
pois.HMM.pn2pw <- function(m, lambda, gamma, delta = NULL,
                           stationary = TRUE) {
  tlambda <- log(lambda)
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if (stationary) {
    # If tdelta is set to NULL and returned in the list,
    # it will cause issues when optimizing with TMB
    return(list(tlambda = tlambda, tgamma = tgamma))
  } else {
    tdelta <- log(delta[- 1] / delta[1])
    # TMB requires a list
    return(list(tlambda = tlambda, tgamma = tgamma, tdelta = tdelta))
  }
}

## ---- pois.HMM.pw2pn
# Transform Poisson working parameters to natural parameters
pois.HMM.pw2pn <- function(m, parvect, stationary = TRUE) {
  parvect <- unlist(parvect)
  lambda <- exp(parvect[1:m])
  gamma <- diag(m)
  if (m == 1) return(list(lambda = lambda, gamma = gamma, delta = 1))
  gamma[!gamma] <- exp(parvect[(m + 1):(m * m)])
  gamma <- gamma / apply(gamma, 1, sum)
  if (stationary) {
    # The code from Zucchini can crash when dealing with computationally small numbers.
    delta <- stat.dist(gamma)
  } else {
    foo <- c(1, exp(parvect[(m * m + 1):(m * m + m - 1)]))
    delta <- foo / sum(foo)
  }
  return(list(lambda = lambda, gamma = gamma, delta = delta))
}

## ---- pois.HMM.generate.sample
# Generate a random sample from a HMM
pois.HMM.generate.sample  <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  }
  x <- rpois(ns, lambda = mod$lambda[state])
  return(list(data = x, state = state))
}

## ---- pois.HMM.generate.estimable.sample
# Generate a random sample from a HMM
pois.HMM.generate.estimable.sample <- function(ns,
                                               mod,
                                               testing_params,
                                               params_names = PARAMS_NAMES,
                                               std_error = FALSE) {
  if(anyNA(c(ns, mod, testing_params))) {
    stop("Some parameters are missing in pois.HMM.generate.estimable.sample")
  }
  # Count occurrences for each error
  failure <- c("state_number" = 0,
               "TMB_null" = 0,
               "TMB_converge" = 0,
               "TMB_G_null" = 0,
               "TMB_G_converge" = 0,
               "TMB_H_null" = 0,
               "TMB_H_converge" = 0,
               "TMG_GH_null" = 0,
               "TMG_GH_converge" = 0,
               "NA_value" = 0)
  m <- mod$m
  # Loop as long as there is an issue with nlminb
  repeat {
    mod_temp <- NULL
    #simulate the data
    new_data <- pois.HMM.generate.sample(ns = ns,
                                         mod = mod)
    
    # If the number of states generated is different from m, discard the data
    if (length(unique(new_data$state)) != m) {
      failure["state_number"] <- failure["state_number"] + 1
      next
    }
    
    TMB_new_data <- list(x = new_data$data, m = m)
    
    testing_w_params <- pois.HMM.pn2pw(m = m,
                                       lambda = testing_params$lambda,
                                       gamma = testing_params$gamma,
                                       delta = testing_params$delta)
    
    # Test TMB
    suppressWarnings(mod_temp <- TMB.estimate(TMB_data = TMB_new_data,
                                              parameters = testing_w_params,
                                              std_error = std_error))
    # If nlminb doesn't reach any result, discard the data
    if (is.null(mod_temp)) {
      failure["TMB_null"] <- failure["TMB_null"] + 1
      next
    }
    # If nlminb doesn't converge successfully, discard the data
    if (mod_temp$convergence != 0) {
      failure["TMB_converge"] <- failure["TMB_converge"] + 1
      next
    }
    
    # Test TMB_G
    suppressWarnings(mod_temp <- TMB.estimate(TMB_data = TMB_new_data,
                                              parameters = testing_w_params,
                                              gradient = TRUE,
                                              std_error = std_error))
    # If nlminb doesn't reach any result, discard the data
    if (is.null(mod_temp)) {
      failure["TMB_G_null"] <- failure["TMB_G_null"] + 1
      next
    }
    # If nlminb doesn't converge successfully, discard the data
    if (mod_temp$convergence != 0) {
      failure["TMB_G_converge"] <- failure["TMB_G_converge"] + 1
      next
    }
    
    # Test TMB_H
    suppressWarnings(mod_temp <- TMB.estimate(TMB_data = TMB_new_data,
                                              parameters = testing_w_params,
                                              hessian = TRUE,
                                              std_error = std_error))
    # If nlminb doesn't reach any result, discard the data
    if (is.null(mod_temp)) {
      failure["TMB_H_null"] <- failure["TMB_H_null"] + 1
      next
    }
    # If nlminb doesn't converge successfully, discard the data
    if (mod_temp$convergence != 0) {
      failure["TMB_H_converge"] <- failure["TMB_H_converge"] + 1
      next
    }
    
    # Test TMB_GH
    suppressWarnings(mod_temp <- TMB.estimate(TMB_data = TMB_new_data,
                                              parameters = testing_w_params,
                                              gradient = TRUE,
                                              hessian = TRUE,
                                              std_error = std_error))
    # If nlminb doesn't reach any result, discard the data
    if (is.null(mod_temp)) {
      failure["TMB_GH_null"] <- failure["TMB_GH_null"] + 1
      next
    }
    # If nlminb doesn't converge successfully, discard the data
    if (mod_temp$convergence != 0) {
      failure["TMB_GH_converge"] <- failure["TMB_GH_converge"] + 1
      next
    }
    
    natural_parameters <- list(m = m,
                               lambda = mod_temp$lambda,
                               gamma = mod_temp$gamma,
                               delta = mod_temp$delta,
                               lambda_std_error = mod_temp$lambda_std_error,
                               gamma_std_error = mod_temp$gamma_std_error,
                               delta_std_error = mod_temp$delta_std_error)
    
    # If some parameters are NA for some reason, discard the data
    if (anyNA(natural_parameters[params_names], recursive = TRUE)) {
      failure["NA_value"] <- failure["NA_value"] + 1
      next
    }
    
    # If everything went well, end the "repeat" loop
    break
  }
  return(list(data = new_data$data,
              states = new_data$state,
              natural_parameters = natural_parameters,
              mod = mod_temp,
              failure = failure))
}

## ---- pois.HMM.label.order
# Relabel states by increasing Poisson means
pois.HMM.label.order <- function(m,
                                 lambda,
                                 gamma,
                                 delta = NULL,
                                 lambda_std_error = NULL,
                                 gamma_std_error = NULL,
                                 delta_std_error = NULL) {
  # gamma_vector_indices is used to calculate the indices of the reordered TPM gamma as
  # a vector for reordering the rows of the complete CI data.frame used for the article.
  gamma_vector_indices <- 1:(m ^ 2)
  gamma_vector_matrix <- matrix(gamma_vector_indices, nrow = m, ncol = m)
  ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)
  
  # Get the indices of the sorted states
  # according to ascending lambda
  # sorted_lambda contains the permutations needed
  ordered_lambda_indices <- order(lambda)
  ordered_lambda <- lambda[ordered_lambda_indices]
  names(ordered_lambda) <- NULL
  # Reorder the TPM according to the switched states
  # in the sorted lambda
  ordered_gamma <- matrix(0, nrow = m, ncol = m)
  for (col in 1:m) {
    new_col <- which(ordered_lambda_indices == col)
    for (row in 1:m) {
      new_row <- which(ordered_lambda_indices == row)
      ordered_gamma[row, col] <- gamma[new_row, new_col]
      
      # Reorder the vector TPM
      ordered_gamma_vector_matrix[row, col] <- gamma_vector_matrix[new_row, new_col]
    }
  }
  # Same for the TPM's standard errors
  ordered_gamma_std_error <- NULL
  if (!is.null(gamma_std_error)) {
    ordered_gamma_std_error <- matrix(0, nrow = m, ncol = m)
    for (col in 1:m) {
      new_col <- which(ordered_lambda_indices == col)
      for (row in 1:m) {
        new_row <- which(ordered_lambda_indices == row)
        ordered_gamma_std_error[row, col] <- gamma_std_error[new_row, new_col]
      }
    }
  }
  
  # Reorder the stationary distribution if it is provided
  # Generate it otherwise
  if (is.null(delta)) {
    ordered_delta <- stat.dist(ordered_gamma)
  } else {
    ordered_delta <- delta[ordered_lambda_indices]
  }
  # Reorder the standard errors
  ordered_lambda_std_error <- lambda_std_error[ordered_lambda_indices]
  ordered_delta_std_error <- delta_std_error[ordered_lambda_indices]
  
  # The vector is assumed filled column-wise instead of row-wise,
  # because column-wise is the default way R handles matrix to vector conversion.
  # Change to row-wise if needed by replacing ordered_gamma_vector_matrix with
  # t(ordered_gamma_vector_matrix), or add byrow=TRUE
  # to "ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)"
  # We don't use it in case there is a bug, but it makes logical sense that it should work
  ordered_gamma_vector_matrix <- as.numeric(ordered_gamma_vector_matrix)
  
  result <- list(lambda = ordered_lambda,
                 gamma = ordered_gamma,
                 delta = ordered_delta,
                 lambda_std_error = ordered_lambda_std_error,
                 gamma_std_error = ordered_gamma_std_error,
                 delta_std_error = ordered_delta_std_error,
                 ordered_lambda_indices = ordered_lambda_indices,
                 ordered_gamma_vector_indices = ordered_gamma_vector_matrix,
                 # delta and lambda are the same size, so they are ordered the same way
                 ordered_delta_indices = ordered_lambda_indices)
  
  # Remove the NULL elements
  result[sapply(result, is.null)] <- NULL
  
  return(result)
}

## ---- pois.HMM.mllk
# Calculate the negative log-likelihood, based on the book by Zucchini
pois.HMM.mllk <- function(parvect, x_alias, m_alias, stationary = TRUE) {
  # The variable names m and x are already used as parameters for the hessian
  # m_alias and x_alias are only replacement names
  m <- m_alias
  x <- x_alias
  n <- length(x)
  pn <- pois.HMM.pw2pn(m, parvect, stationary = stationary)
  emission_probs <- get.emission.probs(x, pn$lambda)
  
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
# Compute the ML estimates without using TMB. Bases on pois.HMM.mle in the book by Zucchini.
DM.estimate <- function(x, m, lambda0, gamma0, delta0 = NULL, stationary = TRUE) {
  parvect0 <- pois.HMM.pn2pw(m = m, lambda = lambda0, gamma = gamma0,
                             delta = delta0, stationary = stationary)
  # nlminb needs a vector, not a list
  parvect0 <- unlist(parvect0)
  
  mod <- nlminb(start = parvect0,
                objective = pois.HMM.mllk,
                x_alias = x,
                m_alias = m,
                stationary = stationary)
  pw <- mod$par
  pn <- pois.HMM.pw2pn(m, as.numeric(pw))
  mllk <- mod$objective
  np <- unlist(parvect0)
  np <- length(np)
  AIC <- 2 * (mllk + np)
  n <- sum(!is.na(x))
  BIC <- 2 * mllk + np * log(n)
  return(list(m = m, lambda = pn$lambda, gamma = pn$gamma, delta = pn$delta,
              convergence = mod$convergence, mllk = mllk, AIC = AIC, BIC = BIC,
              mod = mod))
  
}

# Nice publication-worthy theme for ggplot
## ---- theme
theme_Publication <- function(base_size = 10) {
  library(grid)
  library(ggthemes)
  theme_foundation(base_size = base_size) +
    theme(plot.title = element_text(face = "bold",
                                    size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(size = rel(1)),
          axis.title.y = element_text(angle = 90, vjust = 2),
          axis.title.x = element_text(vjust = - 0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour = "#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key.size= unit(0.5, "cm"),
          legend.spacing = unit(0, "cm"),
          legend.title = element_text(),
          plot.margin = unit(c(5, 5, 5, 5), "mm"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "transparent"),
          strip.text = element_text()
    )
}

