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
  # sigma is a array (p rows, p*p columns, m elements)
  # gamma is a matrix (m rows, m columns)
  # delta is a vector (m columns)
  
  # An array of covariance matrices (sigma) gets transformed according to Zucchini (Section 20.1, page 259)
  # This is the log-Cholesky parametrization
  # tsigma is a matrix (m rows, p(p+1)/2 columns)
  tsigma <- apply(X = sigma,
                  MARGIN = 3,
                  FUN = function(e) {
                    res <- chol(e)
                    diag(res) <- log(diag(res))
                    
                    # The (strictly) lower triangular part is 0, so we can ignore it
                    idx_upper_tri <- upper.tri(res, diag = TRUE)
                    
                    # For simplicity with TMB, we vectorize the matrix row-wise
                    result <- res[idx_upper_tri]
                    return(result)
                  })
  
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if (stationary) {
    # If tdelta is set to NULL and returned in the list,
    # it will cause issues when optimizing with TMB
    return(list(tmu = t(mu), tsigma = tsigma, tgamma = tgamma))
  } else {
    tdelta <- log(delta[- 1] / delta[1])
    # TMB requires a list
    return(list(tmu = t(mu), tsigma = tsigma, tgamma = tgamma, tdelta = tdelta))
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
  # x <- rpois(ns, lambda = mod$lambda[state])
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
                                   delta_std_error = NULL) {
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
  # Reorder the standard errors
  ordered_mu_std_error <- mu_std_error[, ordered_mu_indices]
  ordered_sigma_std_error <- sigma_std_error[, , ordered_mu_indices]
  ordered_delta_std_error <- delta_std_error[ordered_mu_indices]
  
  # The vector is assumed filled column-wise instead of row-wise,
  # because column-wise is the default way R handles matrix to vector conversion.
  # Change to row-wise if needed by replacing ordered_gamma_vector_matrix with
  # t(ordered_gamma_vector_matrix), or add byrow=TRUE
  # to "ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)"
  # We don't use it in case there is a bug, but it makes logical sense that it should work
  ordered_gamma_vector_matrix <- as.numeric(ordered_gamma_vector_matrix)
  
  result <- list(mu = ordered_mu,
                 sigma = ordered_sigma,
                 gamma = ordered_gamma,
                 delta = ordered_delta,
                 mu_std_error = ordered_mu_std_error,
                 sigma_std_error = ordered_sigma_std_error,
                 gamma_std_error = ordered_gamma_std_error,
                 delta_std_error = ordered_delta_std_error,
                 ordered_mu_indices = ordered_mu_indices,
                 ordered_sigma_indices = ordered_mu_indices,
                 ordered_gamma_vector_indices = ordered_gamma_vector_matrix,
                 # delta and mu are the same size, so they are ordered the same way
                 ordered_delta_indices = ordered_mu_indices)
  
  # Remove the NULL elements
  result[sapply(result, is.null)] <- NULL
  
  return(result)
}
