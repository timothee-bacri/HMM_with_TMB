# Utility functions for hmm
# The "## ---- function_name" text lets us define code chunks to import and display in the supplement on GitHub easily
# without needing to copy-paste.

## ---- norm.TMB.estimate
# Estimation using TMB, wrapper of the optimizers, capable of moderate error handling
# Requires either TMB_data & working_parameters, or obj
norm.TMB.estimate <- function(TMB_data = NULL,
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
  if (return_MSE_estimators == TRUE && (is.null(TMB_data$true_mu) ||
                                        is.null(TMB_data$true_sigma) ||
                                        is.null(TMB_data$true_gamma) ||
                                        is.null(TMB_data$true_delta))) {
    message("You have asked to return the MSE_estimators but failed to pass all the true parameter values as argument")
    return()
  }
  
  m <- TMB_data$m
  
  if (return_MSE_estimators == TRUE) {
    dll <- "norm_hmm_MSE_estimators"
  } else if (return_smoothing_probabilities == TRUE) {
    dll <- "norm_hmm_smoothing"
    if (is.null(TMB_data$start_row) || is.null(TMB_data$start_col) || is.null(TMB_data$nb_rows) || is.null(TMB_data$nb_cols)) {
      TMB_data$start_row <- NA
      TMB_data$start_col <- NA
      TMB_data$nb_rows <- NA
      TMB_data$nb_cols <- NA
    }
  } else {
    dll <- "norm_hmm"
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
  
  rows <- rownames(adrep) == "mu"
  mu <- adrep[rows, "Estimate"]
  mu_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL

  rows <- rownames(adrep) == "sigma"
  sigma <- adrep[rows, "Estimate"]
  sigma_std_error <- if (std_error == TRUE) adrep[rows, "Std. Error"] else NULL
  
  rows <- rownames(adrep) == "gamma"
  gamma <- adrep[rows, "Estimate"]
  gamma <- matrix(gamma, ncol = m)
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
  
  # Label switching
  # It is necessary in order to bootstrap correctly
  if (label_switch == TRUE) {
    final_result <- norm.HMM.label.order(m = m,
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

## ---- norm.HMM.decode
# Computes log-forward, log-backward and conditional probabilities
# and decoding based on an optimized TMB::MakeADFun object
norm.HMM.decode <- function(obj) {
  
  # Setup
  # Retrieve the objects at ML value
  adrep <- obj$report(obj$env$last.par.best)
  delta <- adrep$delta
  gamma <- adrep$gamma
  mu <- adrep$mu
  sigma <- adrep$sigma
  emission_probs <- get.norm.emission.probs(obj$env$data$x, mu, sigma)
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
  tgamma_indices <- max(sigma_indices) + 1:m
  delta <- max(tgamma_indices) + 1:m
  
  parvect <- unlist(parvect)
  mu <- parvect[mu_indices]
  sigma <- exp(parvect[sigma_indices])
  gamma <- diag(m)
  if (m == 1) return(list(tmu = mu, sigma = sigma, gamma = gamma, delta = 1))
  gamma[!gamma] <- exp(parvect[tgamma_indices])
  gamma <- gamma / apply(gamma, 1, sum)
  if (stationary) {
    # The code from Zucchini can crash when dealing with computationally small numbers.
    delta <- stat.dist(gamma)
  } else {
    foo <- c(1, exp(parvect[(m * m + 1):(m * m + m - 1)]))
    delta <- foo / sum(foo)
  }
  return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
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


## ---- norm.HMM.generate.estimable.sample
# Generate a random sample from a HMM
norm.HMM.generate.estimable.sample <- function(ns,
                                               mod,
                                               testing_params = mod,
                                               params_names = PARAMS_NAMES_NORM,
                                               return_std_error = FALSE,
                                               return_the_MSE_estimators = FALSE,
                                               # label_switch = TRUE,
                                               control_args = CONTROL_ARGS,
                                               debug_message = "") {
  if(anyNA(c(ns, mod, testing_params))) {
    stop("Some parameters are NA in norm.HMM.generate.estimable.sample")
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
    new_data <- norm.HMM.generate.sample(ns = ns,
                                         mod = mod)
    
    # If the number of states generated is different from m, discard the data
    if (length(unique(new_data$state)) != m) {
      problems["state_number"] <- problems["state_number"] + 1
      next
    }
    
    TMB_new_data <- list(x = new_data$data,
                         m = m)
    
    if (return_the_MSE_estimators == TRUE) {
      TMB_new_data$true_mu <- mod$mu
      TMB_new_data$true_sigma <- mod$sigma
      TMB_new_data$true_gamma <- mod$gamma
      TMB_new_data$true_delta <- mod$delta
    }
    
    testing_w_params <- norm.HMM.pn2pw(m = m,
                                       mu = testing_params$mu,
                                       sigma = testing_params$sigma,
                                       gamma = testing_params$gamma,
                                       delta = testing_params$delta)
    
    # If at least an optimizer shows an issue, regenerate a sample
    problem <- FALSE
    
    # Test BFGS
    result_BFGS <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                     working_parameters = testing_w_params,
                                     std_error = return_std_error,
                                     return_MSE_estimators = return_the_MSE_estimators,
                                     optimizer = "BFGS",
                                     debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_BFGS)) {
      problems["BFGS_error"] <- problems["BFGS_error"] + 1
      problem <- TRUE
    } else if (result_BFGS$convergence == FALSE) {
      problems["BFGS_convergence_failures"] <- problems["BFGS_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test BFGS with gradient
    result_BFGS_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                        working_parameters = testing_w_params,
                                        gradient = TRUE,
                                        std_error = return_std_error,
                                        return_MSE_estimators = return_the_MSE_estimators,
                                        optimizer = "BFGS",
                                        debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_BFGS_gr)) {
      problems["BFGS_gr_error"] <- problems["BFGS_gr_error"] + 1
      problem <- TRUE
    } else if (result_BFGS_gr$convergence == FALSE) {
      problems["BFGS_gr_convergence_failures"] <- problems["BFGS_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test L-BFGS-B
    result_L_BFGS_B <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                         working_parameters = testing_w_params,
                                         std_error = return_std_error,
                                         return_MSE_estimators = return_the_MSE_estimators,
                                         optimizer = "L-BFGS-B",
                                         debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_L_BFGS_B)) {
      problems["L_BFGS_B_error"] <- problems["L_BFGS_B_error"] + 1
      problem <- TRUE
    } else if (result_L_BFGS_B$convergence == FALSE) {
      problems["L_BFGS_B_convergence_failures"] <- problems["L_BFGS_B_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test L-BFGS-B with gradient
    result_L_BFGS_B_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                            working_parameters = testing_w_params,
                                            gradient = TRUE,
                                            std_error = return_std_error,
                                            return_MSE_estimators = return_the_MSE_estimators,
                                            optimizer = "L-BFGS-B",
                                            debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_L_BFGS_B_gr)) {
      problems["L_BFGS_B_gr_error"] <- problems["L_BFGS_B_gr_error"] + 1
      problem <- TRUE
    } else if (result_L_BFGS_B_gr$convergence == FALSE) {
      problems["L_BFGS_B_gr_convergence_failures"] <- problems["L_BFGS_B_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test CG
    result_CG <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                   working_parameters = testing_w_params,
                                   std_error = return_std_error,
                                   return_MSE_estimators = return_the_MSE_estimators,
                                   optimizer = "CG",
                                   debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_CG)) {
      problems["CG_error"] <- problems["CG_error"] + 1
      problem <- TRUE
    } else if (result_CG$convergence == FALSE) {
      problems["CG_convergence_failures"] <- problems["CG_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test CG with gradient
    result_CG_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                      working_parameters = testing_w_params,
                                      gradient = TRUE,
                                      std_error = return_std_error,
                                      return_MSE_estimators = return_the_MSE_estimators,
                                      optimizer = "CG",
                                      debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_CG_gr)) {
      problems["CG_gr_error"] <- problems["CG_gr_error"] + 1
      problem <- TRUE
    } else if (result_CG_gr$convergence == FALSE) {
      problems["CG_gr_convergence_failures"] <- problems["CG_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test Nelder-Mead
    result_Nelder_Mead <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                            working_parameters = testing_w_params,
                                            std_error = return_std_error,
                                            return_MSE_estimators = return_the_MSE_estimators,
                                            optimizer = "Nelder-Mead",
                                            debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_Nelder_Mead)) {
      problems["Nelder_Mead_error"] <- problems["Nelder_Mead_error"] + 1
      problem <- TRUE
    } else if (result_Nelder_Mead$convergence == FALSE) {
      problems["Nelder_Mead_convergence_failures"] <- problems["Nelder_Mead_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test nlm
    result_nlm <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                    working_parameters = testing_w_params,
                                    std_error = return_std_error,
                                    return_MSE_estimators = return_the_MSE_estimators,
                                    optimizer = "nlm",
                                    debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlm)) {
      problems["nlm_error"] <- problems["nlm_error"] + 1
      problem <- TRUE
    } else if (result_nlm$convergence == FALSE) {
      problems["nlm_convergence_failures"] <- problems["nlm_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlm with gradient
    result_nlm_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       gradient = TRUE,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "nlm",
                                       debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlm_gr)) {
      problems["nlm_gr_error"] <- problems["nlm_gr_error"] + 1
      problem <- TRUE
    } else if (result_nlm_gr$convergence == FALSE) {
      problems["nlm_gr_convergence_failures"] <- problems["nlm_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlm with hessian
    result_nlm_he <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       hessian = TRUE,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "nlm",
                                       debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlm_he)) {
      problems["nlm_he_error"] <- problems["nlm_he_error"] + 1
      problem <- TRUE
    } else if (result_nlm_he$convergence == FALSE) {
      problems["nlm_he_convergence_failures"] <- problems["nlm_he_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlm with gradient and hessian
    result_nlm_grhe <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                         working_parameters = testing_w_params,
                                         gradient = TRUE,
                                         hessian = TRUE,
                                         std_error = return_std_error,
                                         return_MSE_estimators = return_the_MSE_estimators,
                                         optimizer = "nlm",
                                         debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlm_grhe)) {
      problems["nlm_grhe_error"] <- problems["nlm_grhe_error"] + 1
      problem <- TRUE
    } else if (result_nlm_grhe$convergence == FALSE) {
      problems["nlm_grhe_convergence_failures"] <- problems["nlm_grhe_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test nlminb
    result_nlminb <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "nlminb",
                                       debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb)) {
      problems["nlminb_error"] <- problems["nlminb_error"] + 1
      problem <- TRUE
    } else if (result_nlminb$convergence == FALSE) {
      problems["nlminb_convergence_failures"] <- problems["nlminb_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlminb with gradient
    result_nlminb_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                          working_parameters = testing_w_params,
                                          gradient = TRUE,
                                          std_error = return_std_error,
                                          return_MSE_estimators = return_the_MSE_estimators,
                                          optimizer = "nlminb",
                                          debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb_gr)) {
      problems["nlminb_gr_error"] <- problems["nlminb_gr_error"] + 1
      problem <- TRUE
    } else if (result_nlminb_gr$convergence == FALSE) {
      problems["nlminb_gr_convergence_failures"] <- problems["nlminb_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlminb with hessian
    result_nlminb_he <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                          working_parameters = testing_w_params,
                                          hessian = TRUE,
                                          std_error = return_std_error,
                                          return_MSE_estimators = return_the_MSE_estimators,
                                          optimizer = "nlminb",
                                          debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb_he)) {
      problems["nlminb_he_error"] <- problems["nlminb_he_error"] + 1
      problem <- TRUE
    } else if (result_nlminb_he$convergence == FALSE) {
      problems["nlminb_he_convergence_failures"] <- problems["nlminb_he_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test nlminb with gradient and hessian
    result_nlminb_grhe <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                            working_parameters = testing_w_params,
                                            gradient = TRUE,
                                            hessian = TRUE,
                                            std_error = return_std_error,
                                            return_MSE_estimators = return_the_MSE_estimators,
                                            optimizer = "nlminb",
                                            debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_nlminb_grhe)) {
      problems["nlminb_grhe_error"] <- problems["nlminb_grhe_error"] + 1
      problem <- TRUE
    } else if (result_nlminb_grhe$convergence == FALSE) {
      problems["nlminb_grhe_convergence_failures"] <- problems["nlminb_grhe_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test hjn
    result_hjn <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                    working_parameters = testing_w_params,
                                    std_error = return_std_error,
                                    return_MSE_estimators = return_the_MSE_estimators,
                                    optimizer = "hjn",
                                    debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_hjn)) {
      problems["hjn_error"] <- problems["hjn_error"] + 1
      problem <- TRUE
    } else if (result_hjn$convergence == FALSE) {
      problems["hjn_convergence_failures"] <- problems["hjn_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test marqLevAlg
    result_marqLevAlg <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                           working_parameters = testing_w_params,
                                           std_error = return_std_error,
                                           return_MSE_estimators = return_the_MSE_estimators,
                                           optimizer = "marqLevAlg",
                                           debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg)) {
      problems["marqLevAlg_error"] <- problems["marqLevAlg_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg$convergence == FALSE) {
      problems["marqLevAlg_convergence_failures"] <- problems["marqLevAlg_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test marqLevAlg with gradient
    result_marqLevAlg_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                              working_parameters = testing_w_params,
                                              gradient = TRUE,
                                              std_error = return_std_error,
                                              return_MSE_estimators = return_the_MSE_estimators,
                                              optimizer = "marqLevAlg",
                                              debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg_gr)) {
      problems["marqLevAlg_gr_error"] <- problems["marqLevAlg_gr_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg_gr$convergence == FALSE) {
      problems["marqLevAlg_gr_convergence_failures"] <- problems["marqLevAlg_gr_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test marqLevAlg with hessian
    result_marqLevAlg_he <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                              working_parameters = testing_w_params,
                                              hessian = TRUE,
                                              std_error = return_std_error,
                                              return_MSE_estimators = return_the_MSE_estimators,
                                              optimizer = "marqLevAlg",
                                              debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg_he)) {
      problems["marqLevAlg_he_error"] <- problems["marqLevAlg_he_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg_he$convergence == FALSE) {
      problems["marqLevAlg_he_convergence_failures"] <- problems["marqLevAlg_he_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # Test marqLevAlg with gradient and hessian
    result_marqLevAlg_grhe <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                                working_parameters = testing_w_params,
                                                gradient = TRUE,
                                                hessian = TRUE,
                                                std_error = return_std_error,
                                                return_MSE_estimators = return_the_MSE_estimators,
                                                optimizer = "marqLevAlg",
                                                debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_marqLevAlg_grhe)) {
      problems["marqLevAlg_grhe_error"] <- problems["marqLevAlg_grhe_error"] + 1
      problem <- TRUE
    } else if (result_marqLevAlg_grhe$convergence == FALSE) {
      problems["marqLevAlg_grhe_convergence_failures"] <- problems["marqLevAlg_grhe_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # # Test ucminf
    # result_ucminf <- norm.TMB.estimate(TMB_data = TMB_new_data,
    #                             working_parameters = testing_w_params,
    #                             gradient = TRUE,
    #                             hessian = TRUE,
    #                             std_error = return_std_error,
    #                             return_MSE_estimators = return_the_MSE_estimators,
    #                             optimizer = "ucminf",
    #                             debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    # if (is.null(result_ucminf)) {
    #   problems["ucminf_error"] <- problems["ucminf_error"] + 1
    #   problem <- TRUE
    # } else if (result_ucminf$convergence == FALSE) {
    #   problems["ucminf_convergence_failures"] <- problems["ucminf_convergence_failures"] + 1
    #   problem <- TRUE
    # }
    
    
    # Test newuoa
    result_newuoa <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                       working_parameters = testing_w_params,
                                       std_error = return_std_error,
                                       return_MSE_estimators = return_the_MSE_estimators,
                                       optimizer = "newuoa",
                                       debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_newuoa)) {
      problems["newuoa_error"] <- problems["newuoa_error"] + 1
      problem <- TRUE
    } else if (result_newuoa$convergence == FALSE) {
      problems["newuoa_convergence_failures"] <- problems["newuoa_convergence_failures"] + 1
      problem <- TRUE
    }
    
    
    # Test BBoptim
    result_BBoptim <- norm.TMB.estimate(TMB_data = TMB_new_data,
                                        working_parameters = testing_w_params,
                                        std_error = return_std_error,
                                        return_MSE_estimators = return_the_MSE_estimators,
                                        optimizer = "BBoptim",
                                        debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    if (is.null(result_BBoptim)) {
      problems["BBoptim_error"] <- problems["BBoptim_error"] + 1
      problem <- TRUE
    } else if (result_BBoptim$convergence == FALSE) {
      problems["BBoptim_convergence_failures"] <- problems["BBoptim_convergence_failures"] + 1
      problem <- TRUE
    }
    
    # # Test BBoptim with gradient
    # result_BBoptim_gr <- norm.TMB.estimate(TMB_data = TMB_new_data,
    #                                   working_parameters = testing_w_params,
    #                                   gradient = TRUE,
    #                                   std_error = return_std_error,
    #                                   return_MSE_estimators = return_the_MSE_estimators,
    #                                   optimizer = "BBoptim",
    #                                   debug_message = paste(debug_message, "> norm.HMM.generate.estimable.sample"))
    # if (is.null(result_BBoptim_gr)) {
    #   problems["BBoptim_gr_error"] <- problems["BBoptim_gr_error"] + 1
    #   problem <- TRUE
    # } else if (result_BBoptim_gr$convergence == FALSE) {
    #   problems["BBoptim_gr_convergence_failures"] <- problems["BBoptim_gr_convergence_failures"] + 1
    #   problem <- TRUE
    # }
    
    # Retrieve nlminb's parameters
    natural_parameters <- list(m = m,
                               mu = result_nlminb$mu,
                               sigma = result_nlminb$sigma,
                               gamma = result_nlminb$gamma,
                               delta = result_nlminb$delta,
                               mu_std_error = result_nlminb$mu_std_error,
                               sigma_std_error = result_nlminb$sigma_std_error,
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
  
  rownames_mles <- c(paste0("mu", 1:m), paste0("sigma", 1:m), paste0("gamma", 1:(m ^ 2)), paste0("delta", 1:m))
  
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

## ---- norm.HMM.label.order
# Relabel states by increasing Gaussian means
norm.HMM.label.order <- function(m,
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
  # according to ascending mu
  # ordered_mu contains the permutations needed
  ordered_mu_indices <- order(mu)
  ordered_mu <- mu[ordered_mu_indices]
  
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
  
  ordered_sigma <- sigma[ordered_mu_indices]
  # Reorder the stationary distribution if it is provided
  # Generate it otherwise
  ordered_delta <- if (!is.null(delta)) delta[ordered_mu_indices] else stat.dist(ordered_gamma)
  
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
  ordered_mu_std_error <- mu_std_error[ordered_mu_indices]
  ordered_sigma_std_error <- sigma_std_error[ordered_mu_indices]
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
  
  # Remove the NULL elements
  result[sapply(result, is.null)] <- NULL
  
  return(result)
}

## ---- norm.HMM.mllk
# Calculate the negative log-likelihood, based on the book by Zucchini
norm.HMM.mllk <- function(parvect, x_alias, m_alias, stationary = TRUE) {
  # The variable names m and x are already used as parameters for the hessian
  # m_alias and x_alias are only replacement names
  m <- m_alias
  x <- x_alias
  n <- length(x)
  pn <- norm.HMM.pw2pn(m = m, parvect = parvect, stationary = stationary)
  emission_probs <- get.norm.emission.probs(data = x, mu = pn$mu, sigma = pn$sigma)
  
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

## ---- norm.latex.names
# Get latex names for parameters (mu, sigma, gamma, delta)
norm.latex.names <- function(m) {
  params_names_latex <- paste0(rep("$\\mu_{",
                                   m),
                               1:m,
                               "}$")
  params_names_latex <- c(params_names_latex,
                          paste0(rep("$\\sigma_{",
                                     m),
                                 1:m,
                                 "}$"))
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
