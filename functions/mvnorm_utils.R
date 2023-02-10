# Utility functions for multivariate Gaussian hmm
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

## ---- mvnorm.get.emission.probs
# Calculate emission probabilities
mvnorm.get.emission.probs <- function(data, mu, sigma) {
  n <- length(data[, 1])
  m <- nrow(mu)
  emission_probs <- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) {
    if (is.na(data[i])) {
      emission_probs[i, ] <- rep(1, m)
    } else {
      for (m_idx in 1:m) {
        emission_probs[i, m_idx] <- mvtnorm::dmvnorm(x = data[i, ],
                                                     mean = mu[m_idx, ],
                                                     sigma = sigma[, , m_idx])
      }
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

## ---- mvnorm.HMM.pn2pw
# Transform multivariate Gaussian natural parameters to working parameters
mvnorm.HMM.pn2pw <- function(m, mu, sigma, gamma, delta = NULL,
                             stationary = TRUE) {
  # mu is a matrix (m rows, p columns) (doesn't require transformation)
  # sigma is an array (p rows, p columns, m matrices)
  # gamma is a matrix (m rows, m columns)
  # delta is a vector (m columns)
  
  # An array of covariance matrices (sigma) gets transformed according to Zucchini (Section 20.1, page 259)
  # This is the log-Cholesky parametrization
  # tsigma is a matrix (p(p+1)/2 rows, m columns)
  tsigma <- apply(X = sigma,
                  MARGIN = 3,
                  FUN = function(e) {
                    res <- chol(e)
                    diag(res) <- log(diag(res))
                    
                    # The (strictly) lower triangular part is 0, so we can ignore it
                    idx_upper_tri <- upper.tri(res, diag = TRUE)
                    
                    # For simplicity with TMB, we vectorize the matrix column-wise
                    result <- res[idx_upper_tri]
                    return(result)
                  })
  
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  # TMB requires a list
  if (stationary) {
    # If tdelta is set to NULL and returned in the list,
    # it will cause issues when optimizing with TMB
    return(list(tmu = mu, tsigma = tsigma, tgamma = tgamma))
  } else {
    tdelta <- log(delta[- 1] / delta[1])
    return(list(tmu = mu, tsigma = tsigma, tgamma = tgamma, tdelta = tdelta))
  }
}

## ---- mvnorm.HMM.pw2pn
# Transform multivariate Gaussian working parameters to natural parameters
mvnorm.HMM.pw2pn <- function(m, p, parvect, stationary = TRUE) {
  tmu_indices <- 1:(m * p)
  tsigma_indices <- max(tmu_indices) + 1:(m * p * (p + 1) / 2)
  tgamma_indices <- max(tsigma_indices) + 1:(m ^ 2 - m)
  delta_indices <- max(tgamma_indices) + 1:m
  
  parvect <- unlist(parvect)
  
  tmu <- parvect[tmu_indices]
  mu <- matrix(tmu, nrow = p, ncol = m)
  
  tsigma <- parvect[tsigma_indices]
  sigma_array <- array(0, dim = c(p, p, m))
  
  # Fill upper triangular elements with working parameters column-wise
  tsigma_p_indices <- 0
  for (m_idx in 1:m) {
    tsigma_p_indices <- max(tsigma_p_indices) + 1:(p * (p + 1) / 2)
    sigma_array_indices <- upper.tri(sigma_array[, , m_idx], diag = TRUE)
    sigma_array[, , m_idx][sigma_array_indices] <- tsigma[tsigma_p_indices]
    
    # Undo the log transformation of the diagonal
    diag(sigma_array[, , m_idx]) <- exp(diag(sigma_array[, , m_idx]))
    
    # Undo the Cholesky transformation
    sigma_matrix <- sigma_array[, , m_idx]
    temporary_matrix <- t(sigma_matrix) %*% sigma_matrix
    sigma_array[, , m_idx] <- temporary_matrix
  }
  
  gamma <- diag(m)
  if (m == 1) return(list(tmu = mu, sigma = sigma, gamma = gamma, delta = 1))
  gamma[!gamma] <- exp(parvect[tgamma_indices])
  gamma <- gamma / apply(gamma, 1, sum)
  if (stationary) {
    # The code from Zucchini can crash when dealing with computationally small numbers.
    delta <- stat.dist(gamma)
  } else {
    foo <- c(1, exp(parvect[delta_indices]))
    delta <- foo / sum(foo)
  }
  return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
}

## ---- mvnorm.HMM.generate.sample
# Generate a random sample from a HMM
mvnorm.HMM.generate.sample  <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  }
  x <- matrix(NA, nrow = ns, ncol = p)
  for (i in 1:ns) {
    x[i, ] <- rmvnorm(n = 1,
                      mean = mod$mu[state[i], ],
                      sigma = mod$sigma[, , state[i]])
    
  }
  names(x) <- NULL
  return(list(data = x, state = state))
}

## ---- mvnorm.HMM.label.order
# Relabel states by first index increasing Gaussian means
mvnorm.HMM.label.order <- function(m,
                                   mu,
                                   sigma,
                                   gamma,
                                   delta = NULL,
                                   mu_std_error = NULL,
                                   sigma_std_error = NULL,
                                   gamma_std_error = NULL,
                                   delta_std_error = NULL,
                                   smoothing_probs = NULL,
                                   smoothing_probs_std_error = NULL,
                                   ldecode = NULL,
                                   indices = FALSE) {
  # gamma_vector_indices is used to calculate the indices of the reordered TPM gamma as
  # a vector for reordering the rows of the complete CI data.frame used for the article.
  gamma_vector_indices <- 1:(m ^ 2)
  gamma_vector_matrix <- matrix(gamma_vector_indices, nrow = m, ncol = m)
  ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)
  
  # Get the indices of the sorted states
  # according to ascending mu[, 1] (first index of the m vectors mu of size p)
  # ordered_mu contains the permutations needed
  ordered_mu_indices <- order(mu[, 1])
  ordered_mu <- mu[, ordered_mu_indices]
  ordered_sigma <- sigma[, , ordered_mu_indices]
  # names(ordered_mu) <- NULL
  
  # Reorder the TPM according to the switched states
  # in the sorted mu
  ordered_gamma <- matrix(0, nrow = m, ncol = m)
  for (col in 1:m) {
    new_col <- which(ordered_mu_indices == col)
    for (row in 1:m) {
      new_row <- which(ordered_mu_indices == row)
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
      new_col <- which(ordered_mu_indices == col)
      for (row in 1:m) {
        new_row <- which(ordered_mu_indices == row)
        ordered_gamma_std_error[row, col] <- gamma_std_error[new_row, new_col]
      }
    }
  }
  
  # Reorder the stationary distribution if it is provided
  # Generate it otherwise
  if (is.null(delta)) {
    ordered_delta <- stat.dist(ordered_gamma)
  } else {
    ordered_delta <- delta[ordered_mu_indices]
  }  
  
  # Reorder the smoothing probabilities (1 row per state, 1 column per data)
  ordered_smoothing_probs <- if (!is.null(smoothing_probs)) smoothing_probs[ordered_mu_indices, ] else NULL
  # Reorder the decoded states
  if (!is.null(ldecode)) {
    ldecode <- sapply(ldecode, function(e) {
      # ldecode = -1 for missing data
      # Keep -1 in that case
      if (is.na(e)) {
        return(NA)
      } else {
        return(which(ordered_mu_indices == e))
      }
    })
  } else {
    ldecode <- NULL
  }
  
  # Reorder the standard errors
  ordered_mu_std_error <- mu_std_error[ordered_mu_indices, ]
  ordered_sigma_std_error <- sigma_std_error[, , ordered_mu_indices]
  ordered_delta_std_error <- delta_std_error[ordered_mu_indices]
  
  # 1 row per state, 1 column per data
  ordered_smoothing_probs_std_error <- if (!is.null(smoothing_probs_std_error)) smoothing_probs_std_error[ordered_mu_indices, ] else NULL
  
  # The vector is assumed filled column-wise instead of row-wise,
  # because column-wise is the default way R handles matrix to vector conversion.
  # Change to row-wise if needed by replacing ordered_gamma_vector_matrix with
  # t(ordered_gamma_vector_matrix), or add byrow=TRUE
  # to "ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)"
  # We don't use it in case there is a bug, but it makes logical sense that it should work
  ordered_gamma_vector_matrix <- as.numeric(ordered_gamma_vector_matrix)
  
  if (indices == TRUE) {
    result <- list(m = m,
                   mu = ordered_mu,
                   gamma = ordered_gamma,
                   delta = ordered_delta,
                   mu_std_error = ordered_mu_std_error,
                   sigma_std_error = ordered_sigma_std_error,
                   gamma_std_error = ordered_gamma_std_error,
                   delta_std_error = ordered_delta_std_error,
                   smoothing_probs = ordered_smoothing_probs,
                   smoothing_probs_std_error = ordered_smoothing_probs_std_error,
                   ldecode = ldecode,
                   ordered_mu_indices = ordered_mu_indices,
                   ordered_gamma_vector_indices = ordered_gamma_vector_matrix,
                   # delta, sigma and mu are the same size, so they are ordered the same way
                   ordered_sigma_indices = ordered_mu_indices,
                   ordered_delta_indices = ordered_mu_indices)
  } else {
    result <- list(m = m,
                   mu = ordered_mu,
                   sigma = ordered_sigma,
                   gamma = ordered_gamma,
                   delta = delta,
                   mu_std_error = ordered_mu_std_error,
                   sigma_std_error = ordered_sigma_std_error,
                   gamma_std_error = ordered_gamma_std_error,
                   delta_std_error = ordered_delta_std_error,
                   smoothing_probs = ordered_smoothing_probs,
                   smoothing_probs_std_error = ordered_smoothing_probs_std_error,
                   ldecode = ldecode)
  }
}


## ---- mvnorm.HMM.label.order.sigma
# Relabel states by increasing or decreasing (depending on 
# sigma_ordering_decreasing) Gaussian covariance matrix value:
# sigma[sigma_ordering_row, sigma_ordering_column]
mvnorm.HMM.label.order.sigma <- function(m,
                                         mu,
                                         sigma,
                                         gamma,
                                         delta = NULL,
                                         sigma_ordering_row = 1,
                                         sigma_ordering_column = 1,
                                         sigma_ordering_decreasing = FALSE,
                                         mu_std_error = NULL,
                                         sigma_std_error = NULL,
                                         gamma_std_error = NULL,
                                         delta_std_error = NULL,
                                         smoothing_probs = NULL,
                                         smoothing_probs_std_error = NULL,
                                         ldecode = NULL,
                                         indices = FALSE) {
  if (anyNA(c(m, mu, sigma, gamma))) {
    return(NA)
  }
  
  # Remove vector names (optional, but looks better without redundancy)
  names(mu) <- NULL
  names(mu_std_error) <- NULL
  names(sigma) <- NULL
  names(sigma_std_error) <- NULL
  names(delta) <- NULL
  names(delta_std_error) <- NULL
  
  # gamma_vector_indices is used to calculate the indices of the reordered TPM gamma as
  # a vector for reordering the rows of the complete CI data.frame used for the article.
  gamma_vector_indices <- 1:(m ^ 2)
  gamma_vector_matrix <- matrix(gamma_vector_indices, nrow = m, ncol = m)
  ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)
  
  # Get the indices of the sorted states
  # according to ascending sigma[sigma_ordering_row, sigma_ordering_column]
  # ordered_sigma contains the correct ordering
  if (is.null(sigma_ordering_row) || is.null(sigma_ordering_column)) {
    stop(paste("You need to specify which row and column indices to use for covariance ordering in mvnorm.HMM.label.order.sigma"))
  }
  if (! ((sigma_ordering_row %in% 1:m) && (sigma_ordering_column %in% 1:m))) {
    stop(paste("The row or column you specified is out of bounds in mvnorm.HMM.label.order.sigma"))
  }
  sigma_values_to_sort <- apply(X = sigma,
                                MARGIN = 3,
                                FUN = function(cov) {
                                  return(cov[sigma_ordering_row, sigma_ordering_column])
                                })
  ordered_sigma_indices <- order(sigma_values_to_sort,
                                 decreasing = sigma_ordering_decreasing)
  ordered_mu <- mu[ordered_sigma_indices, ]
  ordered_sigma <- sigma[, , ordered_sigma_indices]
  # names(ordered_mu) <- NULL
  
  # Reorder the TPM according to the switched states
  # in the sorted mu
  ordered_gamma <- matrix(0, nrow = m, ncol = m)
  for (col in 1:m) {
    new_col <- which(ordered_sigma_indices == col)
    for (row in 1:m) {
      new_row <- which(ordered_sigma_indices == row)
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
      new_col <- which(ordered_sigma_indices == col)
      for (row in 1:m) {
        new_row <- which(ordered_sigma_indices == row)
        ordered_gamma_std_error[row, col] <- gamma_std_error[new_row, new_col]
      }
    }
  }
  
  # Reorder the stationary distribution if it is provided
  # Generate it otherwise
  ordered_delta <- if (!is.null(delta)) delta[ordered_sigma_indices] else stat.dist(ordered_gamma)
  
  # Reorder the smoothing probabilities (1 row per state, 1 column per data)
  ordered_smoothing_probs <- if (!is.null(smoothing_probs)) smoothing_probs[ordered_sigma_indices, ] else NULL
  # Reorder the decoded states
  if (!is.null(ldecode)) {
    ldecode <- sapply(ldecode, function(e) {
      # ldecode = -1 for missing data
      # Keep -1 in that case
      if (is.na(e)) {
        return(NA)
      } else {
        return(which(ordered_sigma_indices == e))
      }
    })
  } else {
    ldecode <- NULL
  }
  
  # Reorder the standard errors
  ordered_mu_std_error <- mu_std_error[ordered_sigma_indices, ]
  ordered_sigma_std_error <- sigma_std_error[, , ordered_sigma_indices]
  ordered_delta_std_error <- delta_std_error[ordered_sigma_indices]
  
  # 1 row per state, 1 column per data
  ordered_smoothing_probs_std_error <- if (!is.null(smoothing_probs_std_error)) smoothing_probs_std_error[ordered_sigma_indices, ] else NULL
  
  # The vector is assumed filled column-wise instead of row-wise,
  # because column-wise is the default way R handles matrix to vector conversion.
  # Change to row-wise if needed by replacing ordered_gamma_vector_matrix with
  # t(ordered_gamma_vector_matrix), or add byrow=TRUE
  # to "ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)"
  # We don't use it in case there is a bug, but it makes logical sense that it should work
  ordered_gamma_vector_matrix <- as.numeric(ordered_gamma_vector_matrix)
  
  if (indices == TRUE) {
    result <- list(m = m,
                   mu = ordered_mu,
                   gamma = ordered_gamma,
                   delta = ordered_delta,
                   mu_std_error = ordered_mu_std_error,
                   sigma_std_error = ordered_sigma_std_error,
                   gamma_std_error = ordered_gamma_std_error,
                   delta_std_error = ordered_delta_std_error,
                   smoothing_probs = ordered_smoothing_probs,
                   smoothing_probs_std_error = ordered_smoothing_probs_std_error,
                   ldecode = ldecode,
                   ordered_mu_indices = ordered_sigma_indices,
                   ordered_gamma_vector_indices = ordered_gamma_vector_matrix,
                   # delta, sigma and mu are the same size, so they are ordered the same way
                   ordered_sigma_indices = ordered_sigma_indices,
                   ordered_delta_indices = ordered_sigma_indices)
  } else {
    result <- list(m = m,
                   mu = ordered_mu,
                   sigma = ordered_sigma,
                   gamma = ordered_gamma,
                   delta = delta,
                   mu_std_error = ordered_mu_std_error,
                   sigma_std_error = ordered_sigma_std_error,
                   gamma_std_error = ordered_gamma_std_error,
                   delta_std_error = ordered_delta_std_error,
                   smoothing_probs = ordered_smoothing_probs,
                   smoothing_probs_std_error = ordered_smoothing_probs_std_error,
                   ldecode = ldecode)
  }
  
  # Remove the NULL elements
  result[sapply(result, is.null)] <- NULL
  
  return(result)
}


## ---- norm.TMB.estimate
# Estimation using TMB, wrapper of the optimizers, capable of moderate error handling
# Requires either TMB_data & working_parameters, or obj
mvnorm.TMB.estimate <- function(TMB_data = NULL,
                                working_parameters = NULL,
                                map = list(),
                                obj = NULL,
                                gradient = FALSE,
                                hessian = FALSE,
                                std_error = FALSE,
                                control_list = list(BFGS = list(trace = 0,
                                                                maxit = 10000),

                                                    "L-BFGS-B" = list(trace = 0,
                                                                      maxit = 10000),
                                                    
                                                    CG = list(trace = 0,
                                                              maxit = 10000),

                                                    "Nelder-Mead" = list(trace = 0,
                                                                         maxit = 10000),

                                                    nlm = list(iterlim = 10000),
                                                    
                                                    nlminb = list(iter.max = 10000),

                                                    marqLevAlg = list(maxiter = 10000),
                                                    
                                                    BBoptim = list(trace = FALSE,
                                                                   maxit = 10000)),
                                return_mod = FALSE,
                                return_obj = FALSE,
                                return_MSE_estimators = FALSE,
                                return_smoothing_probabilities = FALSE,
                                optimizer = "nlminb",
                                supported_optimizers = c("BFGS",
                                                         "L-BFGS-B",
                                                         "CG",
                                                         "Nelder-Mead",
                                                         "nlm",
                                                         "nlminb",
                                                         "hjn",
                                                         "marqLevAlg",
                                                         "newuoa",
                                                         "BBoptim"),
                                label_switch = TRUE,
                                label_switch_order_sigma = FALSE,
                                debug_message = "") {
  
  if (!(optimizer %in% supported_optimizers)) {
    message("TMB optimizer argument invalid")
    return()
  }
  
  m <- TMB_data$m
  p <- ncol(TMB_data$x)
  
  if (return_smoothing_probabilities == TRUE) {
    dll <- "mvnorm_hmm_smoothing_truncated"
    if (is.null(TMB_data$start_row) || is.null(TMB_data$start_col) || is.null(TMB_data$nb_rows) || is.null(TMB_data$nb_cols)) {
      TMB_data$start_row <- NA
      TMB_data$start_col <- NA
      TMB_data$nb_rows <- NA
      TMB_data$nb_cols <- NA
    }
  } else {
    dll <- "mvnorm_hmm"
  }
  obj <- if (!is.null(obj)) obj else MakeADFun(TMB_data, working_parameters, DLL = dll, silent = TRUE, map = map)
  
  gr <- if (gradient) obj$gr else NULL
  he <- if (hessian) obj$he else NULL
  
  # Optimizing while handling errors
  mod <- tryCatch(
    # Default control arguments?
    if (is.null(control_list)) {
      # switch is cleaner than many if else conditions
      switch(
        optimizer,
        "BFGS" = optim(par = obj$par,
                       fn = obj$fn,
                       gr = gr,
                       method = "BFGS"),
        "L-BFGS-B" = optim(par = obj$par,
                           fn = obj$fn,
                           gr = gr,
                           method = "L-BFGS-B"),
        "CG" = optim(par = obj$par,
                     fn = obj$fn,
                     gr = gr,
                     method = "CG"),
        "Nelder-Mead" = optim(par = obj$par,
                              fn = obj$fn,
                              method = "Nelder-Mead"),
        "nlm" = nlm(f = nlm_gradient_hessian_objective,
                    p = obj$par,
                    fn = obj$fn,
                    gr = gr,
                    he = he),
        "nlminb" = nlminb(start = obj$par,
                          objective = obj$fn,
                          gradient = gr,
                          hessian = he),
        "hjn" = hjn(par = obj$par,
                    fn = obj$fn),
        "marqLevAlg" = marqLevAlg(b = obj$par,
                                  fn = obj$fn,
                                  gr = gr,
                                  hess = he),
        # "ucminf" = ucminf(par = obj$par,
        #                         fn = obj$fn,
        #                         gr = obj$gr
        # ),
        "newuoa" = newuoa(par = obj$par,
                          fn = obj$fn),
        # The warnings for this optimizer are not supposed to show up, but they do because of a package coding mistake
        "BBoptim" = suppressWarnings(BBoptim(par = obj$par,
                                             fn = obj$fn,
                                             gr = gr,
                                             quiet = TRUE,
                                             control = list(trace = FALSE)))
      )
    } else {
      switch(
        optimizer,
        "BFGS" = optim(par = obj$par,
                       fn = obj$fn,
                       gr = gr,
                       method = "BFGS",
                       control = control_list$BFGS # default list()
        ),
        "L-BFGS-B" = optim(par = obj$par,
                           fn = obj$fn,
                           gr = gr,
                           method = "L-BFGS-B",
                           control = control_list$`L-BFGS-B` # default list()
        ),
        "CG" = optim(par = obj$par,
                     fn = obj$fn,
                     gr = gr,
                     method = "CG",
                     control = control_list$CG # default list()
        ),
        "Nelder-Mead" = optim(par = obj$par,
                              fn = obj$fn,
                              method = "Nelder-Mead",
                              control = control_list$`Nelder-Mead` # default list()
        ),
        "nlm" = nlm(f = nlm_gradient_hessian_objective,
                    p = obj$par,
                    fn = obj$fn,
                    gr = gr,
                    he = he,
                    iterlim = control_list$nlm$iterlim # default 100
        ),
        "nlminb" = nlminb(start = obj$par,
                          objective = obj$fn,
                          gradient = gr,
                          hessian = he,
                          control = control_list$nlminb # default list()
        ),
        "hjn" = hjn(par = obj$par,
                    fn = obj$fn
                    # lower = control_list$hjn$lower,
                    # upper = control_list$hjn$upper
        ),
        "marqLevAlg" = marqLevAlg(b = obj$par,
                                  fn = obj$fn,
                                  gr = gr,
                                  hess = he,
                                  maxiter = control_list$marqLevAlg$maxiter # default 500
        ),
        # "ucminf" = ucminf(par = obj$par,
        #                         fn = obj$fn,
        #                         gr = obj$gr,
        #                         control = control_list$ucminf
        # ),
        "newuoa" = newuoa(par = obj$par,
                          fn = obj$fn
                          # control = control_list$newuoa
        ),
        # The warnings for this optimizer are not supposed to show up, but they do because of a package coding mistake
        "BBoptim" = suppressWarnings(BBoptim(par = obj$par,
                                             fn = obj$fn,
                                             gr = gr,
                                             quiet = TRUE,
                                             # control = list(trace = FALSE)))
                                             control = control_list$BBoptim) # default list()
        )
      )
    },
    error = function(e) {
      message(optimizer, " error with TMB.estimate:")
      message("m = ", m)
      message("gradient = ", gradient)
      message("hessian = ", hessian)
      message("control_args ", if (is.null(control_list)) "our own" else "optimizer default")
      message("debug_message = ", paste(debug_message, "> TMB.estimate"))
      message("The original error message is:\n", e)
      return(NULL)
    }
  )
  
  if (is.null(mod)) {
    message(optimizer, " returns NULL with TMB.estimate:")
    message("m = ", m)
    message("gradient = ", gradient)
    message("hessian = ", hessian)
    message("debug_message = ", paste(debug_message, "> TMB.estimate"))
    return(NULL)
  }
  
  convergence <- switch(
    optimizer,
    "BFGS" = mod$convergence == 0,
    "BFGS_gr" = mod$convergence == 0,
    "L-BFGS-B" = mod$convergence == 0,
    "L-BFGS-B_gr" = mod$convergence == 0,
    "CG" = mod$convergence == 0,
    "CG_gr" = mod$convergence == 0,
    "Nelder-Mead" = mod$convergence == 0,
    "nlm" = mod$code %in% 1:2,
    "nlm_gr" = mod$code %in% 1:2,
    "nlm_he" = mod$code %in% 1:2,
    "nlm_grhe" = mod$code %in% 1:2,
    "nlminb" = mod$convergence == 0,
    "nlminb_gr" = mod$convergence == 0,
    "nlminb_he" = mod$convergence == 0,
    "nlminb_grhe" = mod$convergence == 0,
    "hjn" = mod$convergence == 0,
    "marqLevAlg" = mod$istop == 1,
    "marqLevAlg_gr" = mod$istop == 1,
    "marqLevAlg_he" = mod$istop == 1,
    "marqLevAlg_grhe" = mod$istop == 1,
    # "ucminf" = mod$convergence %in% c(1, 2, 4),
    "newuoa" = mod$ierr == 0,
    "BBoptim" = mod$convergence == 0
    # "BBoptim_gr" = mod$convergence == 0
  )
  if (convergence == FALSE) {
    w <- paste0(optimizer, " convergence failure with TMB.estimate:\n",
                "The original warning message is:", mod$message,
                "\nm = ", m,
                "\ngradient = ", gradient,
                "\nhessian = ", hessian,
                "\ndebug_message = ", paste(debug_message, "> TMB.estimate"),
                "\n")
    warning(w)
  }
  
  # Retrieve results
  adrep <- summary(sdreport(obj), "report")
  
  rows <- rownames(adrep) == "mllk"
  if (length(rows != 0)) {
    nll <- adrep[rows, "Estimate"]
    nll_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  } else {
    nll <- obj$fn(obj$env$last.par.best)
    nll_std_error <- NULL
    message(paste0("Add \"ADREPORT(mllk);\" to ", dll, ".cpp to retrieve the nll's standard error. ", debug_message))
  }
  
  np <- length(unlist(working_parameters))
  AIC <- 2 * (nll + np)
  n <- sum(!is.na(TMB_data$x))
  BIC <- 2 * nll + np * log(n)
  
  rows <- rownames(adrep) == "mu"
  mu <- adrep[rows, "Estimate"]
  mu <- matrix(mu, ncol = p)
  mu_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  mu_std_error <- if (std_error == TRUE) matrix(mu_std_error, ncol = p) else NULL
  
  rows <- rownames(adrep) == "sigma"
  sigma <- adrep[rows, "Estimate"]
  # There are m matrices of dimension p * p, read column-wise
  sigma_array <- array(dim = c(p, p, m))
  indices_sigma_i <- 1:(p^2)
  for (i in 1:m) {
    sigma_i <- sigma[indices_sigma_i]
    sigma_i <- matrix(sigma_i, nrow = p, ncol = p)
    sigma_array[, , i] <- sigma_i
    
    indices_sigma_i <- indices_sigma_i + p^2
  }
  sigma <- sigma_array
  
  if (std_error == TRUE) {
    sigma_std_error <- adrep[rows, "Std. Error"]
    
    sigma_std_error_array <- array(dim = c(p, p, m))
    indices_sigma_std_error_i <- 1:(p^2)
    for (i in 1:m) {
      sigma_std_error_i <- sigma_std_error[indices_sigma_std_error_i]
      sigma_std_error_i <- matrix(sigma_std_error_i, nrow = p, ncol = p)
      sigma_std_error_array[, , i] <- sigma_std_error_i
      
      indices_sigma_std_error_i <- indices_sigma_std_error_i + p^2
    }
    sigma_std_error <- sigma_std_error_array
  } else {
    sigma_std_error <- NULL
  }
  
  rows <- rownames(adrep) == "gamma"
  gamma <- adrep[rows, "Estimate"]
  gamma <- matrix(gamma, ncol = m)
  gamma_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  gamma_std_error <- if (std_error == TRUE) matrix(gamma_std_error, nrow = m, ncol = m) else NULL
  
  rows <- rownames(adrep) == "delta"
  delta <- adrep[rows, "Estimate"]
  delta_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  
  smoothing_probs <- NULL
  smoothing_probs_std_error <- NULL
  ldecode <- NULL
  if (return_smoothing_probabilities == TRUE) {
    rows <- rownames(adrep) == "truncated_smoothing_probs"
    # We follow the notation of (Zucchini et al., 2016)
    # One row per hidden state, one column per data
    smoothing_probs <- adrep[rows, "Estimate"]
    smoothing_probs <- matrix(smoothing_probs, nrow = m)
    smoothing_probs_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
    smoothing_probs_std_error <- if (std_error == TRUE) matrix(smoothing_probs_std_error, nrow = m) else NULL
    
    ldecode <- obj$report(obj$env$last.par.best)$ldecode
  }
  
  iterations <- switch(
    optimizer,
    "BFGS" = NA,
    "L-BFGS-B" = NA,
    "CG" = NA,
    "Nelder-Mead" = NA,
    "nlm" = mod$iterations,
    "nlminb" = mod$iterations,
    "hjn" = NA,
    "marqLevAlg" = mod$ni,
    # "ucminf" = mod$info["neval"],
    "newuoa" = NA,
    "BBoptim" = mod$iter
  )
  
  # Label switching
  # It is necessary in order to bootstrap correctly
  if (label_switch == TRUE) {
    if (label_switch_order_sigma == TRUE) {
      final_result <- mvnorm.HMM.label.order.sigma(m = m,
                                                   mu = mu,
                                                   sigma = sigma,
                                                   gamma = gamma,
                                                   delta = delta,
                                                   mu_std_error = mu_std_error,
                                                   sigma_std_error = sigma_std_error,
                                                   gamma_std_error = gamma_std_error,
                                                   delta_std_error = delta_std_error,
                                                   smoothing_probs = smoothing_probs,
                                                   smoothing_probs_std_error = smoothing_probs_std_error,
                                                   ldecode = ldecode)
    } else {
      final_result <- mvnorm.HMM.label.order(m = m,
                                             mu = mu,
                                             sigma = sigma,
                                             gamma = gamma,
                                             delta = delta,
                                             mu_std_error = mu_std_error,
                                             sigma_std_error = sigma_std_error,
                                             gamma_std_error = gamma_std_error,
                                             delta_std_error = delta_std_error,
                                             smoothing_probs = smoothing_probs,
                                             smoothing_probs_std_error = smoothing_probs_std_error,
                                             ldecode = ldecode)
    }
  } else {
    final_result <- list(m = m,
                         mu = mu,
                         sigma = sigma,
                         gamma = gamma,
                         delta = delta,
                         mu_std_error = mu_std_error,
                         sigma_std_error = sigma_std_error,
                         gamma_std_error = gamma_std_error,
                         delta_std_error = delta_std_error,
                         # MSE_estimators = MSE_estimators,
                         # MSE_estimators_std_error = MSE_estimators_std_error,
                         smoothing_probs = smoothing_probs,
                         smoothing_probs_std_error = smoothing_probs_std_error,
                         ldecode = ldecode)
  }
  
  final_result["nll"] <- nll
  final_result["nll_std_error"] <- nll_std_error
  final_result["convergence"] <- convergence
  final_result["iterations"] <- iterations
  final_result["AIC"] <- AIC
  final_result["BIC"] <- BIC
  # obj and mod are large objects to display and return. Don't return them unless explicitly asked.
  if (return_obj == FALSE) final_result <- append(final_result, list(mod = mod))
  if (return_mod == FALSE) final_result <- append(final_result, list(obj = obj))
  
  # Remove NULL results
  final_result[sapply(final_result, is.null)] <- NULL
  return(final_result)
}