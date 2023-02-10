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
# Estimation using TMB, wrapper of the optimizers, capable of moderate error handling
# Requires either TMB_data & working_parameters, or obj
TMB.estimate <- function(TMB_data = NULL,
                         working_parameters = NULL,
                         map = list(),
                         obj = NULL,
                         gradient = FALSE,
                         hessian = FALSE,
                         std_error = FALSE,
                         control_list = CONTROL_ARGS,
                         return_mod = FALSE,
                         return_obj = FALSE,
                         return_MSE_estimators = FALSE,
                         return_smoothing_probabilities = FALSE,
                         optimizer = "nlminb",
                         supported_optimizers = OPTIMIZERS_METHOD,
                         label_switch = TRUE,
                         debug_message = "") {
  
  if (!(optimizer %in% supported_optimizers)) {
    message("TMB optimizer argument invalid")
    return()
  }
  if (return_MSE_estimators == TRUE && (is.null(TMB_data$true_lambda) ||
                                        is.null(TMB_data$true_gamma) ||
                                        is.null(TMB_data$true_delta))) {
    message("You have asked to return the MSE_estimators but failed to pass all the true parameter values as argument")
    return()
  }
  
  m <- TMB_data$m
  
  if (return_MSE_estimators == TRUE) {
    dll <- "poi_hmm_MSE_estimators"
  } else if (return_smoothing_probabilities == TRUE) {
    dll <- "poi_hmm_smoothing"
    if (is.null(TMB_data$start_row) || is.null(TMB_data$start_col) || is.null(TMB_data$nb_rows) || is.null(TMB_data$nb_cols)) {
      TMB_data$start_row <- NA
      TMB_data$start_col <- NA
      TMB_data$nb_rows <- NA
      TMB_data$nb_cols <- NA
    }
  } else {
    dll <- "poi_hmm"
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
  nll <- adrep[rows, "Estimate"]
  nll_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  
  np <- length(unlist(working_parameters))
  AIC <- 2 * (nll + np)
  n <- sum(!is.na(TMB_data$x))
  BIC <- 2 * nll + np * log(n)
  
  rows <- rownames(adrep) == "lambda"
  lambda <- adrep[rows, "Estimate"]
  lambda_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  
  rows <- rownames(adrep) == "gamma"
  gamma <- adrep[rows, "Estimate"]
  gamma <- matrix(gamma, nrow = m, ncol = m)
  gamma_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  gamma_std_error <- if (std_error == TRUE) matrix(gamma_std_error, nrow = m, ncol = m) else NULL
  
  rows <- rownames(adrep) == "delta"
  delta <- adrep[rows, "Estimate"]
  delta_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  
  MSE_estimators <- NULL
  MSE_estimators_std_error <- NULL
  if (return_MSE_estimators == TRUE) {
    rows <- rownames(adrep) == "MSE_estimators"
    MSE_estimators <- adrep[rows, "Estimate"]
    MSE_estimators_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  }
  
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
  
  # obj and mod are large objects to display and return. Don't return them unless explicitly asked
  if (return_obj == FALSE) obj <- NULL
  if (return_mod == FALSE) mod <- NULL
  
  final_result <- list(m = m,
                       lambda = lambda,
                       gamma = gamma,
                       delta = delta,
                       lambda_std_error = lambda_std_error,
                       gamma_std_error = gamma_std_error,
                       delta_std_error = delta_std_error,
                       # MSE_estimators = MSE_estimators,
                       # MSE_estimators_std_error = MSE_estimators_std_error,
                       smoothing_probs = smoothing_probs,
                       smoothing_probs_std_error = smoothing_probs_std_error,
                       ldecode = ldecode)
  # Label switching
  # It is necessary in order to bootstrap correctly
  if (label_switch == TRUE) {
    final_result <- pois.HMM.label.order(m = m,
                                         lambda = lambda,
                                         gamma = gamma,
                                         delta = delta,
                                         lambda_std_error = lambda_std_error,
                                         gamma_std_error = gamma_std_error,
                                         delta_std_error = delta_std_error,
                                         smoothing_probs = smoothing_probs,
                                         smoothing_probs_std_error = smoothing_probs_std_error,
                                         ldecode = ldecode)
  }
  
  final_result["MSE_estimators"] <- MSE_estimators
  final_result["MSE_estimators_std_error"] <- MSE_estimators_std_error
  final_result["nll"] <- nll
  final_result["nll_std_error"] <- nll_std_error
  final_result["convergence"] <- convergence
  final_result["iterations"] <- iterations
  final_result["AIC"] <- AIC
  final_result["BIC"] <- BIC
  final_result <- append(final_result, list(mod = mod))
  final_result <- append(final_result, list(obj = obj))
  
  # Remove NULL results
  final_result[sapply(final_result, is.null)] <- NULL
  return(final_result)
}

## ---- pois.HMM.decode
# Computes log-forward, log-backward and conditional probabilities
# and decoding based on an optimized TMB::MakeADFun object
pois.HMM.decode <- function(obj) {
  
  # Setup
  # Retrieve the objects at ML value
  adrep <- obj$report(obj$env$last.par.best)
  delta <- adrep$delta
  gamma <- adrep$gamma
  lambda <- adrep$lambda
  emission_probs <- get.emission.probs(obj$env$data$x, lambda)
  n <- adrep$n
  m <- length(delta)
  mllk <- obj$fn(obj$env$last.par.best)
  
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

# ## ---- pois.HMM.forecast
# # Given x and mod, equation (5.5) in (Zucchini et al., 2016) is used here to find forecast probabilities
# # h is the number of steps in the future we want to predict
# # The range of x-values for which these probabilities are required is specified by the input xf
# # and the range of times for which they are required by h
# # The model estimates are provided by mod (list containing lambda, gamma, delta)
# # The forecast probabilities are conditional on the past data provided by x (vector)
# pois.HMM.forecast <- function(xf, h = 1, x, mod) {
#   n <- length(x)
#   nxf <- length(xf)
#   dxf <- matrix (0, nrow = h, ncol = nxf)
#   foo <- mod$delta * dpois(x[1], mod$lambda)
#   sumfoo <- sum(foo)
#   lscale <- log(sumfoo)
#   foo <- foo/sumfoo
#   for (i in 2:n)
#   {
#     foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
#     sumfoo <- sum(foo)
#     lscale <- lscale + log(sumfoo)
#     foo <- foo / sumfoo
#   }
#   for (i in 1:h) {
#     foo <- foo %*% mod$gamma
#     for (j in 1:mod$m) dxf[i, ] <- dxf[i, ] + foo[j] * dpois(xf, mod$lambda[j])
#   }
#   return(dxf)
# }

## ---- quantile.colwise
# Computes the 2.5% and 97.5% quantiles
quantile.colwise <- function(data, ...) {
  return(quantile(data, probs = c(0.05 / 2, 1 - 0.05 / 2), ...))
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
                                               testing_params = mod,
                                               params_names = PARAMS_NAMES,
                                               return_std_error = FALSE,
                                               return_the_MSE_estimators = FALSE,
                                               control_args = CONTROL_ARGS,
                                               debug_message = "") {
  if(anyNA(c(ns, mod, testing_params))) {
    stop("Some parameters are NA in pois.HMM.generate.estimable.sample")
  }
  # Count occurrences for each error
  problems <- c("state_number" = 0,
                "BFGS_null" = 0,
                "BFGS_convergence_failures" = 0,
                "BFGS_gr_null" = 0,
                "BFGS_gr_convergence_failures" = 0,
                "L_BFGS_B_null" = 0,
                "L_BFGS_B_convergence_failures" = 0,
                "L_BFGS_B_gr_null" = 0,
                "L_BFGS_B_gr_convergence_failures" = 0,
                "CG_null" = 0,
                "CG_convergence_failures" = 0,
                "CG_gr_null" = 0,
                "CG_gr_convergence_failures" = 0,
                "Nelder_Mead_null" = 0,
                "Nelder_Mead_convergence_failures" = 0,
                "nlm_null" = 0,
                "nlm_convergence_failures" = 0,
                "nlm_gr_null" = 0,
                "nlm_gr_convergence_failures" = 0,
                "nlm_he_null" = 0,
                "nlm_he_convergence_failures" = 0,
                "nlm_grhe_null" = 0,
                "nlm_grhe_convergence_failures" = 0,
                "nlminb_null" = 0,
                "nlminb_convergence_failures" = 0,
                "nlminb_gr_null" = 0,
                "nlminb_gr_convergence_failures" = 0,
                "nlminb_he_null" = 0,
                "nlminb_he_convergence_failures" = 0,
                "nlminb_grhe_null" = 0,
                "nlminb_grhe_convergence_failures" = 0,
                "hjn_null" = 0,
                "hjn_convergence_failures" = 0,
                "marqLevAlg_null" = 0,
                "marqLevAlg_convergence_failures" = 0,
                "marqLevAlg_gr_null" = 0,
                "marqLevAlg_gr_convergence_failures" = 0,
                "marqLevAlg_he_null" = 0,
                "marqLevAlg_he_convergence_failures" = 0,
                "marqLevAlg_grhe_null" = 0,
                "marqLevAlg_grhe_convergence_failures" = 0,
                # "ucminf_null" = 0,
                # "ucminf_convergence_failures" = 0,
                "newuoa_null" = 0,
                "newuoa_convergence_failures" = 0,
                "BBoptim_null" = 0,
                "BBoptim_convergence_failures" = 0,
                # "BBoptim_gr_null" = 0,
                # "BBoptim_gr_convergence_failures" = 0,
                "NA_value" = 0)
  m <- mod$m
  # Loop as long as there is an issue with nlminb, and prevent messages
  # sink(nullfile())
  # sink("log.txt")
  repeat {
    #simulate the data
    new_data <- pois.HMM.generate.sample(ns = ns,
                                         mod = mod)
    
    # If the number of states generated is different from m, discard the data
    if (length(unique(new_data$state)) != m) {
      problems["state_number"] <- problems["state_number"] + 1
      next
    }
    
    TMB_new_data <- list(x = new_data$data,
                         m = m)
    
    if (return_the_MSE_estimators == TRUE) {
      TMB_new_data$true_lambda <- mod$lambda
      TMB_new_data$true_gamma <- mod$gamma
      TMB_new_data$true_delta <- mod$delta
    }
    
    testing_w_params <- pois.HMM.pn2pw(m = m,
                                       lambda = testing_params$lambda,
                                       gamma = testing_params$gamma,
                                       delta = testing_params$delta)
    
    # If at least an optimizer shows an issue, regenerate a sample
    problem <- FALSE
    
    # Test BFGS
    result_BFGS <- TMB.estimate(TMB_data = TMB_new_data,
                                working_parameters = testing_w_params,
                                std_error = return_std_error,
                                return_MSE_estimators = return_the_MSE_estimators,
                                optimizer = "BFGS",
                                debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_BFGS)) {
      problems["BFGS_error"] <- problems["BFGS_error"] + 1
      problem <- TRUE
    } else if (result_BFGS$convergence == FALSE) {
      problems["BFGS_convergence_failures"] <- problems["BFGS_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test BFGS with gradient
    result_BFGS_gr <- TMB.estimate(TMB_data = TMB_new_data,
                                   working_parameters = testing_w_params,
                                   gradient = TRUE,
                                   std_error = return_std_error,
                                   return_MSE_estimators = return_the_MSE_estimators,
                                   optimizer = "BFGS",
                                   debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_BFGS_gr)) {
      problems["BFGS_gr_error"] <- problems["BFGS_gr_error"] + 1
      problem <- TRUE
    } else if (result_BFGS_gr$convergence == FALSE) {
      problems["BFGS_gr_convergence_failures"] <- problems["BFGS_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test L-BFGS-B
    result_L_BFGS_B <- TMB.estimate(TMB_data = TMB_new_data,
                                    working_parameters = testing_w_params,
                                    std_error = return_std_error,
                                    return_MSE_estimators = return_the_MSE_estimators,
                                    optimizer = "L-BFGS-B",
                                    debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_L_BFGS_B)) {
      problems["L_BFGS_B_error"] <- problems["L_BFGS_B_error"] + 1
      problem <- TRUE
    } else if (result_L_BFGS_B$convergence == FALSE) {
      problems["L_BFGS_B_convergence_failures"] <- problems["L_BFGS_B_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test L-BFGS-B with gradient
    result_L_BFGS_B_gr <- TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       gradient = TRUE,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "L-BFGS-B",
                                       debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_L_BFGS_B_gr)) {
      problems["L_BFGS_B_gr_error"] <- problems["L_BFGS_B_gr_error"] + 1
      problem <- TRUE
    } else if (result_L_BFGS_B_gr$convergence == FALSE) {
      problems["L_BFGS_B_gr_convergence_failures"] <- problems["L_BFGS_B_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test CG
    result_CG <- TMB.estimate(TMB_data = TMB_new_data,
                              working_parameters = testing_w_params,
                              std_error = return_std_error,
                              return_MSE_estimators = return_the_MSE_estimators,
                              optimizer = "CG",
                              debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_CG)) {
      problems["CG_error"] <- problems["CG_error"] + 1
      problem <- TRUE
    } else if (result_CG$convergence == FALSE) {
      problems["CG_convergence_failures"] <- problems["CG_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test CG with gradient
    result_CG_gr <- TMB.estimate(TMB_data = TMB_new_data,
                                 working_parameters = testing_w_params,
                                 gradient = TRUE,
                                 std_error = return_std_error,
                                 return_MSE_estimators = return_the_MSE_estimators,
                                 optimizer = "CG",
                                 debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_CG_gr)) {
      problems["CG_gr_error"] <- problems["CG_gr_error"] + 1
      problem <- TRUE
    } else if (result_CG_gr$convergence == FALSE) {
      problems["CG_gr_convergence_failures"] <- problems["CG_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test Nelder-Mead
    result_Nelder_Mead <- TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "Nelder-Mead",
                                       debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_Nelder_Mead)) {
      problems["Nelder_Mead_error"] <- problems["Nelder_Mead_error"] + 1
      problem <- TRUE
    } else if (result_Nelder_Mead$convergence == FALSE) {
      problems["Nelder_Mead_convergence_failures"] <- problems["Nelder_Mead_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test nlm
    result_nlm <- TMB.estimate(TMB_data = TMB_new_data,
                               working_parameters = testing_w_params,
                               std_error = return_std_error,
                               return_MSE_estimators = return_the_MSE_estimators,
                               optimizer = "nlm",
                               debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlm)) {
      problems["nlm_error"] <- problems["nlm_error"] + 1
      problem <- TRUE
    } else if (result_nlm$convergence == FALSE) {
      problems["nlm_convergence_failures"] <- problems["nlm_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlm with gradient
    result_nlm_gr <- TMB.estimate(TMB_data = TMB_new_data,
                                  working_parameters = testing_w_params,
                                  gradient = TRUE,
                                  std_error = return_std_error,
                                  return_MSE_estimators = return_the_MSE_estimators,
                                  optimizer = "nlm",
                                  debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlm_gr)) {
      problems["nlm_gr_error"] <- problems["nlm_gr_error"] + 1
      problem <- TRUE
    } else if (result_nlm_gr$convergence == FALSE) {
      problems["nlm_gr_convergence_failures"] <- problems["nlm_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlm with hessian
    result_nlm_he <- TMB.estimate(TMB_data = TMB_new_data,
                                  working_parameters = testing_w_params,
                                  hessian = TRUE,
                                  std_error = return_std_error,
                                  return_MSE_estimators = return_the_MSE_estimators,
                                  optimizer = "nlm",
                                  debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlm_he)) {
      problems["nlm_he_error"] <- problems["nlm_he_error"] + 1
      problem <- TRUE
    } else if (result_nlm_he$convergence == FALSE) {
      problems["nlm_he_convergence_failures"] <- problems["nlm_he_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlm with gradient and hessian
    result_nlm_grhe <- TMB.estimate(TMB_data = TMB_new_data,
                                    working_parameters = testing_w_params,
                                    gradient = TRUE,
                                    hessian = TRUE,
                                    std_error = return_std_error,
                                    return_MSE_estimators = return_the_MSE_estimators,
                                    optimizer = "nlm",
                                    debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlm_grhe)) {
      problems["nlm_grhe_error"] <- problems["nlm_grhe_error"] + 1
      problem <- TRUE
    } else if (result_nlm_grhe$convergence == FALSE) {
      problems["nlm_grhe_convergence_failures"] <- problems["nlm_grhe_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test nlminb
    result_nlminb <- TMB.estimate(TMB_data = TMB_new_data,
                                  working_parameters = testing_w_params,
                                  std_error = return_std_error,
                                  return_MSE_estimators = return_the_MSE_estimators,
                                  optimizer = "nlminb",
                                  debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb)) {
      problems["nlminb_error"] <- problems["nlminb_error"] + 1
      problem <- TRUE
    } else if (result_nlminb$convergence == FALSE) {
      problems["nlminb_convergence_failures"] <- problems["nlminb_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlminb with gradient
    result_nlminb_gr <- TMB.estimate(TMB_data = TMB_new_data,
                                     working_parameters = testing_w_params,
                                     gradient = TRUE,
                                     std_error = return_std_error,
                                     return_MSE_estimators = return_the_MSE_estimators,
                                     optimizer = "nlminb",
                                     debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb_gr)) {
      problems["nlminb_gr_error"] <- problems["nlminb_gr_error"] + 1
      problem <- TRUE
    } else if (result_nlminb_gr$convergence == FALSE) {
      problems["nlminb_gr_convergence_failures"] <- problems["nlminb_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlminb with hessian
    result_nlminb_he <- TMB.estimate(TMB_data = TMB_new_data,
                                     working_parameters = testing_w_params,
                                     hessian = TRUE,
                                     std_error = return_std_error,
                                     return_MSE_estimators = return_the_MSE_estimators,
                                     optimizer = "nlminb",
                                     debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb_he)) {
      problems["nlminb_he_error"] <- problems["nlminb_he_error"] + 1
      problem <- TRUE
    } else if (result_nlminb_he$convergence == FALSE) {
      problems["nlminb_he_convergence_failures"] <- problems["nlminb_he_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlminb with gradient and hessian
    result_nlminb_grhe <- TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       gradient = TRUE,
                                       hessian = TRUE,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "nlminb",
                                       debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb_grhe)) {
      problems["nlminb_grhe_error"] <- problems["nlminb_grhe_error"] + 1
      problem <- TRUE
    } else if (result_nlminb_grhe$convergence == FALSE) {
      problems["nlminb_grhe_convergence_failures"] <- problems["nlminb_grhe_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test hjn
    result_hjn <- TMB.estimate(TMB_data = TMB_new_data,
                               working_parameters = testing_w_params,
                               std_error = return_std_error,
                               return_MSE_estimators = return_the_MSE_estimators,
                               optimizer = "hjn",
                               debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_hjn)) {
      problems["hjn_error"] <- problems["hjn_error"] + 1
      problem <- TRUE
    } else if (result_hjn$convergence == FALSE) {
      problems["hjn_convergence_failures"] <- problems["hjn_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test marqLevAlg
    result_marqLevAlg <- TMB.estimate(TMB_data = TMB_new_data,
                                      working_parameters = testing_w_params,
                                      std_error = return_std_error,
                                      return_MSE_estimators = return_the_MSE_estimators,
                                      optimizer = "marqLevAlg",
                                      debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg)) {
      problems["marqLevAlg_error"] <- problems["marqLevAlg_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg$convergence == FALSE) {
      problems["marqLevAlg_convergence_failures"] <- problems["marqLevAlg_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test marqLevAlg with gradient
    result_marqLevAlg_gr <- TMB.estimate(TMB_data = TMB_new_data,
                                         working_parameters = testing_w_params,
                                         gradient = TRUE,
                                         std_error = return_std_error,
                                         return_MSE_estimators = return_the_MSE_estimators,
                                         optimizer = "marqLevAlg",
                                         debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg_gr)) {
      problems["marqLevAlg_gr_error"] <- problems["marqLevAlg_gr_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg_gr$convergence == FALSE) {
      problems["marqLevAlg_gr_convergence_failures"] <- problems["marqLevAlg_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test marqLevAlg with hessian
    result_marqLevAlg_he <- TMB.estimate(TMB_data = TMB_new_data,
                                         working_parameters = testing_w_params,
                                         hessian = TRUE,
                                         std_error = return_std_error,
                                         return_MSE_estimators = return_the_MSE_estimators,
                                         optimizer = "marqLevAlg",
                                         debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg_he)) {
      problems["marqLevAlg_he_error"] <- problems["marqLevAlg_he_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg_he$convergence == FALSE) {
      problems["marqLevAlg_he_convergence_failures"] <- problems["marqLevAlg_he_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test marqLevAlg with gradient and hessian
    result_marqLevAlg_grhe <- TMB.estimate(TMB_data = TMB_new_data,
                                           working_parameters = testing_w_params,
                                           gradient = TRUE,
                                           hessian = TRUE,
                                           std_error = return_std_error,
                                           return_MSE_estimators = return_the_MSE_estimators,
                                           optimizer = "marqLevAlg",
                                           debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg_grhe)) {
      problems["marqLevAlg_grhe_error"] <- problems["marqLevAlg_grhe_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg_grhe$convergence == FALSE) {
      problems["marqLevAlg_grhe_convergence_failures"] <- problems["marqLevAlg_grhe_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # # Test ucminf
    # result_ucminf <- TMB.estimate(TMB_data = TMB_new_data,
    #                             working_parameters = testing_w_params,
    #                             gradient = TRUE,
    #                             hessian = TRUE,
    #                             std_error = return_std_error,
    #                             return_MSE_estimators = return_the_MSE_estimators,
    #                             optimizer = "ucminf",
    #                             debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    # if (is.null(result_ucminf)) {
    #   problems["ucminf_error"] <- problems["ucminf_error"] + 1
    #   problem <- TRUE
    # } else if (result_ucminf$convergence == FALSE) {
    #   problems["ucminf_convergence_failures"] <- problems["ucminf_convergence_failures"] + 1
    #   problem <- TRUE
    # }
    
    
    # Test newuoa
    result_newuoa <- TMB.estimate(TMB_data = TMB_new_data,
                                  working_parameters = testing_w_params,
                                  std_error = return_std_error,
                                  return_MSE_estimators = return_the_MSE_estimators,
                                  optimizer = "newuoa",
                                  debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_newuoa)) {
      problems["newuoa_error"] <- problems["newuoa_error"] + 1
      problem <- TRUE
    } else if (result_newuoa$convergence == FALSE) {
      problems["newuoa_convergence_failures"] <- problems["newuoa_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test BBoptim
    result_BBoptim <- TMB.estimate(TMB_data = TMB_new_data,
                                   working_parameters = testing_w_params,
                                   std_error = return_std_error,
                                   return_MSE_estimators = return_the_MSE_estimators,
                                   optimizer = "BBoptim",
                                   debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    if (is.null(result_BBoptim)) {
      problems["BBoptim_error"] <- problems["BBoptim_error"] + 1
      problem <- TRUE
    } else if (result_BBoptim$convergence == FALSE) {
      problems["BBoptim_convergence_failures"] <- problems["BBoptim_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # # Test BBoptim with gradient
    # result_BBoptim_gr <- TMB.estimate(TMB_data = TMB_new_data,
    #                                   working_parameters = testing_w_params,
    #                                   gradient = TRUE,
    #                                   std_error = return_std_error,
    #                                   return_MSE_estimators = return_the_MSE_estimators,
    #                                   optimizer = "BBoptim",
    #                                   debug_message = paste(debug_message, "> pois.HMM.generate.estimable.sample"))
    # if (is.null(result_BBoptim_gr)) {
    #   problems["BBoptim_gr_error"] <- problems["BBoptim_gr_error"] + 1
    #   problem <- TRUE
    # } else if (result_BBoptim_gr$convergence == FALSE) {
    #   problems["BBoptim_gr_convergence_failures"] <- problems["BBoptim_gr_convergence_failures"] + 1
    #   problem <- TRUE
    # }
    
    # Retrieve nlminb's parameters
    natural_parameters <- list(m = m,
                               lambda = result_nlminb$lambda,
                               gamma = result_nlminb$gamma,
                               delta = result_nlminb$delta,
                               lambda_std_error = result_nlminb$lambda_std_error,
                               gamma_std_error = result_nlminb$gamma_std_error,
                               delta_std_error = result_nlminb$delta_std_error)
    
    # If some parameters are NA for some reason, discard the data
    if (anyNA(natural_parameters[params_names], recursive = TRUE)) {
      problems["NA_value"] <- problems["NA_value"] + 1
      problem <- TRUE
    }
    
    # If everything went well, end the "repeat" loop
    # Otherwise, the loop continues
    if (problem == FALSE) {
      break
    }
  }
  
  #sink()
  
  # Retrieve negative log-likelihoods across all optimizers
  all_nlls <- data.frame("BFGS" = result_BFGS$nll,
                         "BFGS_gr" = result_BFGS_gr$nll,
                         "L_BFGS_B" = result_L_BFGS_B$nll,
                         "L_BFGS_B_gr" = result_L_BFGS_B_gr$nll,
                         "CG" = result_CG$nll,
                         "CG_gr" = result_CG_gr$nll,
                         "Nelder_Mead" = result_Nelder_Mead$nll,
                         "nlm" = result_nlm$nll,
                         "nlm_gr" = result_nlm_gr$nll,
                         "nlm_he" = result_nlm_he$nll,
                         "nlm_grhe" = result_nlm_grhe$nll,
                         "nlminb" = result_nlminb$nll,
                         "nlminb_gr" = result_nlminb_gr$nll,
                         "nlminb_he" = result_nlminb_he$nll,
                         "nlminb_grhe" = result_nlminb_grhe$nll,
                         "hjn" = result_hjn$nll,
                         "marqLevAlg" = result_marqLevAlg$nll,
                         "marqLevAlg_gr" = result_marqLevAlg_gr$nll,
                         "marqLevAlg_he" = result_marqLevAlg_he$nll,
                         "marqLevAlg_grhe" = result_marqLevAlg_grhe$nll,
                         # "ucminf" = result_ucminf$nll,
                         "newuoa" = result_newuoa$nll,
                         "BBoptim" = result_BBoptim$nll)
  # "BBoptim_gr" = result_BBoptim_gr$nll)
  
  rownames_mles <- c(paste0("lambda", 1:m), paste0("gamma", 1:(m ^ 2)), paste0("delta", 1:m))
  
  if (return_the_MSE_estimators == TRUE) {
    all_MSE_estimators <- data.frame("BFGS" = result_BFGS$MSE_estimators,
                                     "BFGS_gr" = result_BFGS_gr$MSE_estimators,
                                     "L_BFGS_B" = result_L_BFGS_B$MSE_estimators,
                                     "L_BFGS_B_gr" = result_L_BFGS_B_gr$MSE_estimators,
                                     "CG" = result_CG$MSE_estimators,
                                     "CG_gr" = result_CG_gr$MSE_estimators,
                                     "Nelder_Mead" = result_Nelder_Mead$MSE_estimators,
                                     "nlm" = result_nlm$MSE_estimators,
                                     "nlm_gr" = result_nlm_gr$MSE_estimators,
                                     "nlm_he" = result_nlm_he$MSE_estimators,
                                     "nlm_grhe" = result_nlm_grhe$MSE_estimators,
                                     "nlminb" = result_nlminb$MSE_estimators,
                                     "nlminb_gr" = result_nlminb_gr$MSE_estimators,
                                     "nlminb_he" = result_nlminb_he$MSE_estimators,
                                     "nlminb_grhe" = result_nlminb_grhe$MSE_estimators,
                                     "hjn" = result_hjn$MSE_estimators,
                                     "marqLevAlg" = result_marqLevAlg$MSE_estimators,
                                     "marqLevAlg_gr" = result_marqLevAlg_gr$MSE_estimators,
                                     "marqLevAlg_he" = result_marqLevAlg_he$MSE_estimators,
                                     "marqLevAlg_grhe" = result_marqLevAlg_grhe$MSE_estimators,
                                     # "ucminf" = result_ucminf$MSE_estimators,
                                     "newuoa" = result_newuoa$MSE_estimators,
                                     "BBoptim" = result_BBoptim$MSE_estimators)
    # "BBoptim_gr" = result_BBoptim_gr$MSE_estimators)
  }
  
  all_natural_parameters <- data.frame("m" = m,
                                       "Parameter" = rownames_mles,
                                       "BFGS" = unlist(result_BFGS[params_names]),
                                       "BFGS_gr" = unlist(result_BFGS_gr[params_names]),
                                       "L_BFGS_B" = unlist(result_L_BFGS_B[params_names]),
                                       "L_BFGS_B_gr" = unlist(result_L_BFGS_B_gr[params_names]),
                                       "CG" = unlist(result_CG[params_names]),
                                       "CG_gr" = unlist(result_CG_gr[params_names]),
                                       "Nelder_Mead" = unlist(result_Nelder_Mead[params_names]),
                                       "nlm" = unlist(result_nlm[params_names]),
                                       "nlm_gr" = unlist(result_nlm_gr[params_names]),
                                       "nlm_he" = unlist(result_nlm_he[params_names]),
                                       "nlm_grhe" = unlist(result_nlm_grhe[params_names]),
                                       "nlminb" = unlist(result_nlminb[params_names]),
                                       "nlminb_gr" = unlist(result_nlminb_gr[params_names]),
                                       "nlminb_he" = unlist(result_nlminb_he[params_names]),
                                       "nlminb_grhe" = unlist(result_nlminb_grhe[params_names]),
                                       "hjn" = unlist(result_hjn[params_names]),
                                       "marqLevAlg" = unlist(result_marqLevAlg[params_names]),
                                       "marqLevAlg_gr" = unlist(result_marqLevAlg_gr[params_names]),
                                       "marqLevAlg_he" = unlist(result_marqLevAlg_he[params_names]),
                                       "marqLevAlg_grhe" = unlist(result_marqLevAlg_grhe[params_names]),
                                       # "ucminf" = unlist(result_ucminf[params_names]),
                                       "newuoa" = unlist(result_newuoa[params_names]),
                                       "BBoptim" = unlist(result_BBoptim[params_names]))
  # "BBoptim_gr" = unlist(result_BBoptim_gr[params_names]))
  row.names(all_natural_parameters) <- NULL
  all_natural_parameters_list <- list("BFGS" = result_BFGS[params_names],
                                      "BFGS_gr" = result_BFGS_gr[params_names],
                                      "L_BFGS_B" = result_L_BFGS_B[params_names],
                                      "L_BFGS_B_gr" = result_L_BFGS_B_gr[params_names],
                                      "CG" = result_CG[params_names],
                                      "CG_gr" = result_CG_gr[params_names],
                                      "Nelder_Mead" = result_Nelder_Mead[params_names],
                                      "nlm" = result_nlm[params_names],
                                      "nlm_gr" = result_nlm_gr[params_names],
                                      "nlm_he" = result_nlm_he[params_names],
                                      "nlm_grhe" = result_nlm_grhe[params_names],
                                      "nlminb" = result_nlminb[params_names],
                                      "nlminb_gr" = result_nlminb_gr[params_names],
                                      "nlminb_he" = result_nlminb_he[params_names],
                                      "nlminb_grhe" = result_nlminb_grhe[params_names],
                                      "hjn" = result_hjn[params_names],
                                      "marqLevAlg" = result_marqLevAlg[params_names],
                                      "marqLevAlg_gr" = result_marqLevAlg_gr[params_names],
                                      "marqLevAlg_he" = result_marqLevAlg_he[params_names],
                                      "marqLevAlg_grhe" = result_marqLevAlg_grhe[params_names],
                                      # "ucminf" = result_ucminf[params_names],
                                      "newuoa" = result_newuoa[params_names],
                                      "BBoptim" = result_BBoptim[params_names])
  # "BBoptim_gr" = result_BBoptim_gr[params_names])
  
  all_iterations <- list("BFGS" = result_BFGS$iterations,
                         "BFGS_gr" = result_BFGS_gr$iterations,
                         "L_BFGS_B" = result_L_BFGS_B$iterations,
                         "L_BFGS_B_gr" = result_L_BFGS_B_gr$iterations,
                         "CG" = result_CG$iterations,
                         "CG_gr" = result_CG_gr$iterations,
                         "Nelder_Mead" = result_Nelder_Mead$iterations,
                         "nlm" = result_nlm$iterations,
                         "nlm_gr" = result_nlm_gr$iterations,
                         "nlm_he" = result_nlm_he$iterations,
                         "nlm_grhe" = result_nlm_grhe$iterations,
                         "nlminb" = result_nlminb$iterations,
                         "nlminb_gr" = result_nlminb_gr$iterations,
                         "nlminb_he" = result_nlminb_he$iterations,
                         "nlminb_grhe" = result_nlminb_grhe$iterations,
                         "hjn" = result_hjn$iterations,
                         "marqLevAlg" = result_marqLevAlg$iterations,
                         "marqLevAlg_gr" = result_marqLevAlg_gr$iterations,
                         "marqLevAlg_he" = result_marqLevAlg_he$iterations,
                         "marqLevAlg_grhe" = result_marqLevAlg_grhe$iterations,
                         # "ucminf" = iterations,
                         "newuoa" = result_newuoa$iterations,
                         "BBoptim" = result_BBoptim$iterations)
  
  if (return_std_error == TRUE) {
    params_names <- c(params_names, "nll")
    rownames_mles <- c(rownames_mles, "nll")
    params_names <- paste0(params_names, "_std_error")
    rownames_mles <- paste0(rownames_mles, "_std_error")
    all_std_errors <- data.frame("m" = m,
                                 "Parameter" = rownames_mles,
                                 "BFGS" = unlist(result_BFGS[params_names]),
                                 "BFGS_gr" = unlist(result_BFGS_gr[params_names]),
                                 "L_BFGS_B" = unlist(result_L_BFGS_B[params_names]),
                                 "L_BFGS_B_gr" = unlist(result_L_BFGS_B_gr[params_names]),
                                 "CG" = unlist(result_CG[params_names]),
                                 "CG_gr" = unlist(result_CG_gr[params_names]),
                                 "Nelder_Mead" = unlist(result_Nelder_Mead[params_names]),
                                 "nlm" = unlist(result_nlm[params_names]),
                                 "nlm_gr" = unlist(result_nlm_gr[params_names]),
                                 "nlm_he" = unlist(result_nlm_he[params_names]),
                                 "nlm_grhe" = unlist(result_nlm_grhe[params_names]),
                                 "nlminb" = unlist(result_nlminb[params_names]),
                                 "nlminb_gr" = unlist(result_nlminb_gr[params_names]),
                                 "nlminb_he" = unlist(result_nlminb_he[params_names]),
                                 "nlminb_grhe" = unlist(result_nlminb_grhe[params_names]),
                                 "hjn" = unlist(result_hjn[params_names]),
                                 "marqLevAlg" = unlist(result_marqLevAlg[params_names]),
                                 "marqLevAlg_gr" = unlist(result_marqLevAlg_gr[params_names]),
                                 "marqLevAlg_he" = unlist(result_marqLevAlg_he[params_names]),
                                 "marqLevAlg_grhe" = unlist(result_marqLevAlg_grhe[params_names]),
                                 # "ucminf" = unlist(result_ucminf[params_names]),
                                 "newuoa" = unlist(result_newuoa[params_names]),
                                 "BBoptim" = unlist(result_BBoptim[params_names]))
    # "BBoptim_gr" = unlist(result_BBoptim_gr[params_names]))
    row.names(all_std_errors) <- NULL
  }
  
  result <- list(data = new_data$data,
                 states = new_data$state,
                 natural_parameters = natural_parameters,
                 all_natural_parameters = all_natural_parameters,
                 all_natural_parameters_list = all_natural_parameters_list,
                 all_nlls = all_nlls,
                 all_iterations = all_iterations,
                 problems = problems)
  
  # Return all MSEs derived from estimates
  if (return_the_MSE_estimators == TRUE) {
    result$all_MSE_estimators <- all_MSE_estimators
  }
  
  # Return all std errors
  if (return_std_error == TRUE) {
    result$all_std_errors <- all_std_errors
  }
  
  return(result)
}

## ---- pois.HMM.label.order
# Relabel states by increasing Poisson means
pois.HMM.label.order <- function(m,
                                 lambda,
                                 gamma,
                                 delta = NULL,
                                 lambda_std_error = NULL,
                                 gamma_std_error = NULL,
                                 delta_std_error = NULL,
                                 smoothing_probs = NULL,
                                 smoothing_probs_std_error = NULL,
                                 ldecode = NULL,
                                 indices = FALSE) {
  if (anyNA(c(m, lambda, gamma))) {
    return(NA)
  }
  
  # Remove vector names (optional, but looks better without redundancy)
  names(lambda) <- NULL
  names(lambda_std_error) <- NULL
  names(delta) <- NULL
  names(delta_std_error) <- NULL
  
  # gamma_vector_indices is used to calculate the indices of the reordered TPM gamma as
  # a vector for reordering the rows of the complete CI data.frame used for the article.
  gamma_vector_indices <- 1:(m ^ 2)
  gamma_vector_matrix <- matrix(gamma_vector_indices, nrow = m, ncol = m)
  ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)
  
  # Get the indices of the sorted states
  # according to ascending lambda
  # ordered_lambda contains the permutations needed
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
  ordered_delta <- if (!is.null(delta)) delta[ordered_lambda_indices] else stat.dist(ordered_gamma)
  
  # Reorder the smoothing probabilities (1 row per state, 1 column per data)
  ordered_smoothing_probs <- if (!is.null(smoothing_probs)) smoothing_probs[ordered_lambda_indices, ] else NULL
  # Reorder the decoded states
  if (!is.null(ldecode)) {
    ldecode <- sapply(ldecode, function(e) {
      # ldecode = -1 for missing data
      # Keep -1 in that case
      if (is.na(e)) {
        return(NA)
      } else {
        return(which(ordered_lambda_indices == e))
      }
    })
  } else {
    ldecode <- NULL
  }
  
  # Reorder the standard errors
  ordered_lambda_std_error <- lambda_std_error[ordered_lambda_indices]
  ordered_delta_std_error <- delta_std_error[ordered_lambda_indices]
  # 1 row per state, 1 column per data
  ordered_smoothing_probs_std_error <- if (!is.null(smoothing_probs_std_error)) smoothing_probs_std_error[ordered_lambda_indices, ] else NULL
  
  # The vector is assumed filled column-wise instead of row-wise,
  # because column-wise is the default way R handles matrix to vector conversion.
  # Change to row-wise if needed by replacing ordered_gamma_vector_matrix with
  # t(ordered_gamma_vector_matrix), or add byrow=TRUE
  # to "ordered_gamma_vector_matrix <- matrix(0, nrow = m, ncol = m)"
  # We don't use it in case there is a bug, but it makes logical sense that it should work
  ordered_gamma_vector_matrix <- as.numeric(ordered_gamma_vector_matrix)
  
  if (indices == TRUE) {
    result <- list(m = m,
                   lambda = ordered_lambda,
                   gamma = ordered_gamma,
                   delta = ordered_delta,
                   lambda_std_error = ordered_lambda_std_error,
                   gamma_std_error = ordered_gamma_std_error,
                   delta_std_error = ordered_delta_std_error,
                   smoothing_probs = ordered_smoothing_probs,
                   smoothing_probs_std_error = ordered_smoothing_probs_std_error,
                   ldecode = ldecode,
                   ordered_lambda_indices = ordered_lambda_indices,
                   ordered_gamma_vector_indices = ordered_gamma_vector_matrix,
                   # delta and lambda are the same size, so they are ordered the same way
                   ordered_delta_indices = ordered_lambda_indices)
  } else {
    result <- list(m = m,
                   lambda = ordered_lambda,
                   gamma = ordered_gamma,
                   delta = ordered_delta,
                   lambda_std_error = ordered_lambda_std_error,
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

# ## ---- nlmfn
# # Function to use TMB's gradient and/or hessian in nlm
# nlmfn <- function(par, obj, gr = TRUE, he = TRUE) {
#   res <- as.numeric(obj$fn(par))
#   if(gr) {
#     attr(res, "gradient") <- obj$gr(par)
#   }
#   if(he) {
#     attr(res, "hessian") <- obj$he(par)
#   }
#   return(res)
# }

## ---- nlm_gradient_hessian_objective
# Function to use TMB's gradient and hessian in nlm
nlm_gradient_hessian_objective <- function(par, fn, gr = NULL, he = NULL) {
  res <- fn(par)
  if(!is.null(gr)) attr(res, "gradient") <- gr(par)
  if(!is.null(he)) attr(res, "hessian") <- he(par)
  return(res)
}

## ---- theme
# Nice publication-worthy theme
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
          axis.title = element_text(size = rel(1.2)),
          axis.title.y = element_text(angle = 90, vjust = 2),
          axis.title.x = element_text(vjust = - 0.2),
          # axis.text = element_text(size = rel(1.1)),
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
          # plot.margin = unit(c(5, 5, 5, 5), "mm"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "transparent"),
          strip.text = element_text()
    )
}

## ---- head.secure
# head function that accepts n = 0
head.secure <- function(x, n) {
  if (n == 0) {
    return(x)
  }
  return(head(x, n))
}

## ---- cat.sink
# Safer than sink(nullfile()) ... sink()
# Redirect all console output to oblivion
# closeAllConnections() closes all connections (including sink connections)
sink.safe <- function(...) {
  sink(file = nullfile())
  tryCatch({
    eval(...)
    sink()
  },
  error = function(e) {
    sink()
    print(e)
  })
}

## ---- combine_foreach_rbind
# https://stackoverflow.com/questions/28348749/outptut-two-objects-using-foreach/28354056#28354056
# Combine n data.frames using a list of n data.frames
# Necessary in our case of foreach for parallel programming
combine_foreach_rbind <- function(...) {
  mapply('rbind', ..., SIMPLIFY = FALSE)
}

## ---- pois.HMM.conditional.pseudo.residuals.probs
# Bases on Zucchini (2016) p.82
# Conditional distribution of one data on all other data:
# likelihood of data with u-th observation replaced
# divided by likelihood of data with u-th observation missing
# Observation u is replaced successively by all possibilities in replacement_data_range
pois.HMM.conditional.fitted.probabilities <- function(data, u, replacement_data_range, natural_estimates, stationary = TRUE) {
  # Define a vector, cd, into which the required probabilities will be returned.
  nr <- length(replacement_data_range)
  cd <- rep(NA, nr)
  
  # Transform the parameters into a working format
  m <- length(natural_estimates$lambda)
  w_params <- pois.HMM.pn2pw(m = m, lambda = natural_estimates$lambda, gamma = natural_estimates$gamma, delta = natural_estimates$delta, stationary = stationary)
  
  # Define a vector new_data, which is equal to data but with the u-th
  # observation replaced by NA.
  new_data <- data
  new_data[u] <- NA
  # Use pois.HMM.mllk to compute -log(denominator), i.e. minus one times
  # the log-likelihood of the 'observations', new_data.
  nllk_denomenator <- pois.HMM.mllk(parvect = w_params, m_alias = m, x_alias = new_data, stationary = stationary)
  # Loop over x in replacement_data_range
  for (i in 1:nr) {
    # Now replace the u-th observation by x
    new_data[u] <- replacement_data_range[i]
    # Use pois.HMM.mllk to compute -log(numerator)= -log(likelihood(new_data))
    nllk_numerator <- pois.HMM.mllk(parvect = w_params, m_alias = m, x_alias = new_data, stationary = stationary)
    # Compute the ratio numerator/denominator.
    cd[i] <- exp(- nllk_numerator + nllk_denomenator)
  }
  return(cd)
}

## ---- pois.HMM.fitted.data
# Bases on Zucchini (2016) p.82
# Conditional distribution of one data on all other data:
# likelihood of data with u-th observation replaced
# divided by likelihood of data with u-th observation missing
pois.HMM.conditional.fitted.data <- function(data, u_range, replacement_data_range, natural_estimates, stationary = TRUE) {
  if (max(u_range) > length(data)) {
    stop("The index of u_range in pois.HMM.conditional.fitted.data is too large")
  }
  fitted_data_vector <- c()
  for (u in u_range) {
    probs <- pois.HMM.conditional.fitted.probabilities(data = data, u = u, replacement_data_range = replacement_data_range, natural_estimates = natural_estimates, stationary = stationary)
    max_prob_idx <- which.max(probs)
    fitted_data_point <- replacement_data_range[max_prob_idx]
    
    fitted_data_vector <- c(fitted_data_vector, fitted_data_point)
  }
  return(fitted_data_vector)
}

# source("functions/mclapply.Windows.hack.R")
trunc_custom <- function(x, digits = 0, ...) {
  # base::trunc(x * 10 ^ digits, ...) / 10 ^ digits
  round(x, digits)
}

plot_limits <- function(x) {
  maxi <- max(x, na.rm = TRUE)
  mini <- min(x, na.rm = TRUE)
  if (maxi > 50) {
    return(c(mini, maxi))
  } else if (mini < 0) {
    return(c(mini, maxi + 1))
  } else {
    return(c(0, max(x, 1)))
  }
}

# # Number of decimals on ggplot axis ticks -> 1
# # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2/38722547#38722547
# # How to define the format -> https://www.tutorialspoint.com/c_standard_library/c_function_sprintf.htm
# plot_nb_decimals <- function(x) {
#   maxi <- max(x, na.rm = TRUE)
#   if (maxi > 1000) {
#     return(sprintf("%.0e", x))
#   } else if (50 < maxi && maxi <= 1000) {
#     return(sprintf("%.0f", x))
#   } else if (1 < maxi && maxi <= 50) {
#     return(sprintf("%.1f", x))
#   } else if (maxi <= 1) {
#     return(sprintf("%.2f", x))
#   }
# }

matrix.to.LaTeX <- function(mat) {
  rows <- apply(mat, MARGIN=1, paste, collapse = " & ")
  matrix_string <- paste(rows, collapse = " \\\\ ")
  return(paste("\\begin{pmatrix}", matrix_string, "\\end{pmatrix}"))
}


## ---- pois.latex.names
# Get latex names for parameters (lambda, gamma, delta)
pois.latex.names <- function(m) {
  params_names_latex <- paste0(rep("$\\lambda_{",
                                   m),
                               1:m,
                               "}$")
  for (col in 1:m) {
    # Get row and column indices for gamma instead of the default
    # columnwise index: the default indices are 1:m for the 1st column,
    # then (m + 1):(2 * m) for the 2nd, etc...
    params_names_latex <- c(params_names_latex,
                            paste0(sapply(X = 1:m,
                                          FUN = function(row) {paste0("$\\gamma_{",
                                                                      row,
                                                                      col,
                                                                      "}$")})))
  }
  params_names_latex <- c(params_names_latex,
                          paste0(rep("$\\delta_{",
                                     m),
                                 1:m,
                                 "}$"))
  return(params_names_latex)
}
