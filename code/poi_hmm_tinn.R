FIRST_SEED <- 10
# Prepare the data & parameters, then estimate for different numbers of hidden states
set.seed(FIRST_SEED)
# Parameters and covariates --------------------------
m <- M_LIST_TINN
gamma_init <- matrix(0.2 / (m - 1),
                     nrow = m,
                     ncol = m)
diag(gamma_init) <- 0.8
lambda_init <- seq(quantile(tinn_data,
                            0.1),
                   quantile(tinn_data,
                            0.9),
                   length.out = m)
delta_init <- stat.dist(gamma_init)

# Parameters & covariates for TMB ------------------
working_params_init <- pois.HMM.pn2pw(m,
                                      lambda_init,
                                      gamma_init)
TMB_data <- list(x = tinn_data,
                 m = m)
obj_init <- MakeADFun(TMB_data,
                      working_params_init,
                      DLL = "poi_hmm",
                      silent = TRUE)
parvect_init <- pois.HMM.pn2pw(m = m,
                               lambda = lambda_init,
                               gamma = gamma_init,
                               delta = delta_init)
parvect_init <- unlist(parvect_init)

# Estimation ------------------------------------
dm <- DM.estimate(x = tinn_data,
                  m = m,
                  lambda0 = lambda_init,
                  gamma0 = gamma_init)
tmb <- TMB.estimate(TMB_data = TMB_data,
                    parameters = working_params_init)
tmb_g <- TMB.estimate(TMB_data = TMB_data,
                      parameters = working_params_init,
                      gradient = TRUE)
tmb_h <- TMB.estimate(TMB_data = TMB_data,
                      parameters = working_params_init,
                      hessian = TRUE)
tmb_gh <- TMB.estimate(TMB_data = TMB_data,
                       parameters = working_params_init,
                       gradient = TRUE,
                       hessian = TRUE,
                       std_error = TRUE)

# If one doesn't converge successfully, stop
if (dm$convergence != 0) {
  stop(paste("dm didn't converge properly, tinn dataset, m =",
             m))
}
if (tmb$convergence != 0) {
  stop(paste("tmb didn't converge properly, tinn dataset, m =",
             m))
}
if (tmb_g$convergence != 0) {
  stop(paste("tmb_g didn't converge properly, tinn dataset, m =",
             m))
}
if (tmb_h$convergence != 0) {
  stop(paste("tmb_h didn't converge properly, tinn dataset, m =",
             m))
}
if (tmb_gh$convergence != 0) {
  stop(paste("tmb_gh didn't converge properly, tinn dataset, m =",
             m))
}

# Creating variables for the CIs -----------------
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
len_par <- length(params_names_latex)
indices <- (length(conf_int_tinn$m) + 1):(length(conf_int_tinn$m) + len_par)
conf_int_tinn[indices, "m"] <- m
conf_int_tinn[indices, "Parameter"] <- params_names_latex
# Reminder, PARAMS_NAMES contains c("lambda", "gamma", "delta")
conf_int_tinn[indices, "Estimate"] <- unlist(tmb_gh[PARAMS_NAMES])

param_tmb_gh <- pois.HMM.pn2pw(m = m,
                               lambda = tmb_gh$lambda,
                               gamma = tmb_gh$gamma)

if (m == 1) {
  w_params_names <- c("tlambda1")
} else {
  w_params_names <- c(paste0(rep("tlambda",
                                 m),
                             1:m),
                      paste0(rep("tgamma",
                                 m ^ 2 - m),
                             1:(m ^ 2 - m)))
}
len_w_par <- length(w_params_names)
working_conf_int <- data.frame(w_parameter = w_params_names,
                               lower = rep(NA,
                                           len_w_par),
                               upper = rep(NA,
                                           len_w_par),
                               stringsAsFactors = FALSE)

lambda_indices <- 1:m
gamma_indices <- m + 1:(m ^ 2)
delta_indices <- m ^ 2 + m + (1:m)
tgamma_indices <- (m + 1):(m ^ 2)

# Benchmark the same dataset many times to check the benchmark durations have low variance ---------
set.seed(FIRST_SEED + 1)
if (CONSISTENCY_BENCHMARK_TINN != 0) {
  parvect_benchmark_TMB <- pois.HMM.pn2pw(m = m,
                                          lambda = lambda_init,
                                          gamma = gamma_init,
                                          delta = delta_init)
  TMB_benchmark_data <- TMB_data
  obj_benchmark <- MakeADFun(TMB_benchmark_data,
                             parvect_benchmark_TMB,
                             DLL = "poi_hmm",
                             silent = TRUE)
  # nlminb needs a vector, not a list
  parvect_benchmark_DM <- unlist(parvect_benchmark_TMB)
  
  
  # Parameters & covariates for DM and TMB are already done
  # Estimation benchmark
  temp <- microbenchmark("DM" = nlminb(parvect_benchmark_DM,
                                       pois.HMM.mllk,
                                       x_alias = tinn_data,
                                       m_alias = m)$convergence==0,
                         "TMB" = nlminb(obj_benchmark$par,
                                        obj_benchmark$fn)$convergence==0,
                         "TMB_G" = nlminb(obj_benchmark$par,
                                          obj_benchmark$fn,
                                          gradient = obj_benchmark$gr)$convergence==0,
                         "TMB_H" = nlminb(obj_benchmark$par,
                                          obj_benchmark$fn,
                                          hessian = obj_benchmark$he)$convergence==0,
                         "TMB_GH" = nlminb(obj_benchmark$par,
                                           obj_benchmark$fn,
                                           gradient = obj_benchmark$gr,
                                           hessian = obj_benchmark$he)$convergence==0,
                         times = CONSISTENCY_BENCHMARK_TINN,
                         setup = obj_benchmark <<- MakeADFun(TMB_data,
                                                             parvect_benchmark_TMB,
                                                             DLL = "poi_hmm",
                                                             silent = TRUE),
                         check = "equal")
  
  iterations <- c(nlminb(parvect_benchmark_DM,
                         pois.HMM.mllk,
                         x_alias = tinn_data,
                         m_alias = m)$iterations,
                  nlminb(obj_benchmark$par,
                         obj_benchmark$fn)$iteration,
                  nlminb(obj_benchmark$par,
                         obj_benchmark$fn,
                         gradient = obj_benchmark$gr)$iteration,
                  nlminb(obj_benchmark$par,
                         obj_benchmark$fn,
                         hessian = obj_benchmark$he)$iteration,
                  nlminb(obj_benchmark$par,
                         obj_benchmark$fn,
                         gradient = obj_benchmark$gr,
                         hessian = obj_benchmark$he)$iteration)
  
  times <- temp$time / 10^9
  timeDM <- times[temp$expr == "DM"]
  timeTMB <- times[temp$expr == "TMB"]
  timeTMB_G <- times[temp$expr == "TMB_G"]
  timeTMB_H <- times[temp$expr == "TMB_H"]
  timeTMB_GH <- times[temp$expr == "TMB_GH"]
  consistency_estim_benchmarks_df_tinn <- rbind(
    consistency_estim_benchmarks_df_tinn,
    data.frame(
      time = c(timeDM,
               timeTMB,
               timeTMB_G,
               timeTMB_H,
               timeTMB_GH),
      m = rep(m,
              length(PROCEDURES) * CONSISTENCY_BENCHMARK_TINN),
      procedure = rep(PROCEDURES,
                      each = CONSISTENCY_BENCHMARK_TINN),
      iterations = rep(iterations,
                       each = CONSISTENCY_BENCHMARK_TINN))
  )
}

# Benchmarks ------------
set.seed(FIRST_SEED + 2)
if (BENCHMARK_SAMPLES != 0) {
  for (idx_counter in 1:BENCHMARK_SAMPLES) {
    # Generate data that can be estimated by TMB_GH
    # and is tested on the slightly off parameters from the beginning of this file
    # The goal is to have a dataset that poses no estimation problem
    benchmark_model <- pois.HMM.generate.estimable.sample(
      ns = DATA_SIZE_TINN,
      mod = list(m = m,
                 lambda = tmb_gh$lambda,
                 gamma = tmb_gh$gamma,
                 delta = tmb_gh$delta),
      testing_params = list(m = m,
                            lambda = tmb_gh$lambda,
                            gamma = tmb_gh$gamma,
                            delta = tmb_gh$delta)
    )
    benchmark_data <- benchmark_model$data
    # Benchmark all different combinations of gradient and hessians with DM ----------------
    # Parameters & covariates for DM and TMB
    TMB_benchmark_data <- list(x = benchmark_data,
                               m = m)
    parvect_benchmark_TMB <- pois.HMM.pn2pw(m = m,
                                            lambda = tmb_gh$lambda,
                                            gamma = tmb_gh$gamma,
                                            delta = tmb_gh$delta)
    obj_benchmark <- MakeADFun(TMB_benchmark_data,
                               parvect_benchmark_TMB,
                               DLL = "poi_hmm",
                               silent = TRUE)
    # nlminb needs a vector, not a list
    parvect_benchmark_DM <- unlist(parvect_benchmark_TMB)
    # Estimation benchmark
    temp <- microbenchmark("DM" = nlminb(parvect_benchmark_DM,
                                         pois.HMM.mllk,
                                         x_alias = benchmark_data,
                                         m_alias = m)$convergence==0,
                           "TMB" = nlminb(obj_benchmark$par,
                                          obj_benchmark$fn)$convergence==0,
                           "TMB_G" = nlminb(obj_benchmark$par,
                                            obj_benchmark$fn,
                                            gradient = obj_benchmark$gr)$convergence==0,
                           "TMB_H" = nlminb(obj_benchmark$par,
                                            obj_benchmark$fn,
                                            hessian = obj_benchmark$he)$convergence==0,
                           "TMB_GH" = nlminb(obj_benchmark$par,
                                             obj_benchmark$fn,
                                             gradient = obj_benchmark$gr,
                                             hessian = obj_benchmark$he)$convergence==0,
                           times = 1,
                           check = "equal",
                           setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data,
                                                               parvect_benchmark_TMB,
                                                               DLL = "poi_hmm",
                                                               silent = TRUE))
    times <- temp$time / 10^9
    timeDM <- times[temp$expr == "DM"]
    timeTMB <- times[temp$expr == "TMB"]
    timeTMB_G <- times[temp$expr == "TMB_G"]
    timeTMB_H <- times[temp$expr == "TMB_H"]
    timeTMB_GH <- times[temp$expr == "TMB_GH"]
    
    iterDM <- nlminb(parvect_benchmark_DM,
                     pois.HMM.mllk,
                     x_alias = benchmark_data,
                     m_alias = m)$iterations
    iterTMB <- nlminb(obj_benchmark$par,
                      obj_benchmark$fn)$iterations
    iterTMB_G <- nlminb(obj_benchmark$par,
                        obj_benchmark$fn,
                        gradient = obj_benchmark$gr)$iterations
    iterTMB_H <- nlminb(obj_benchmark$par,
                        obj_benchmark$fn,
                        hessian = obj_benchmark$he)$iterations
    iterTMB_GH <- nlminb(obj_benchmark$par,
                         obj_benchmark$fn,
                         gradient = obj_benchmark$gr,
                         hessian = obj_benchmark$he)$iterations
    
    estim_benchmarks_df_tinn <- rbind(estim_benchmarks_df_tinn,
                                      data.frame(time = c(timeDM,
                                                          timeTMB,
                                                          timeTMB_G,
                                                          timeTMB_H,
                                                          timeTMB_GH),
                                                 m = rep(m,
                                                         length(PROCEDURES)),
                                                 procedure = PROCEDURES,
                                                 iterations = c(iterDM,
                                                                iterTMB,
                                                                iterTMB_G,
                                                                iterTMB_H,
                                                                iterTMB_GH),
                                                 dataset_number = rep(idx_counter,
                                                                      length(PROCEDURES))))
    
  }
}

# Profiling the likelihood --------------------------
for (idx_param in 1:len_w_par) {
  profile <- tmbprofile(obj = tmb_gh$obj,
                        name = idx_param,
                        trace = FALSE)
  
  ci <- tryCatch({
    confint(profile)
  },
  error = function(e){
    return(rep(NA,
               2))
  })
  
  w_param_name <- w_params_names[idx_param]
  working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
  working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
}

# Transform the working parameters into natural ones
# Lambda (m values)
conf_int_tinn$Profile.L[lambda_indices] <- exp(working_conf_int$lower[lambda_indices])
# Gamma (m^2-m working parameters, m^2 natural ones)
if (!anyNA(working_conf_int$lower[tgamma_indices])) {
  natural_gamma <- as.numeric(gamma.w2n(m,
                                        working_conf_int$lower[tgamma_indices]))
  conf_int_tinn$Profile.L[gamma_indices] <- natural_gamma
}
# Lambda (m values)
conf_int_tinn$Profile.U[lambda_indices] <- exp(working_conf_int$upper[lambda_indices])
# Gamma (m^2-m working parameters, m^2 natural ones)
if (!anyNA(working_conf_int$upper[tgamma_indices])) {
  natural_gamma <- as.numeric(gamma.w2n(m,
                                        working_conf_int$upper[tgamma_indices]))
  conf_int_tinn$Profile.U[gamma_indices] <- natural_gamma
}

# Bootstrap ---------------------------
set.seed(FIRST_SEED + 3)
bootstrap_tinn <- data.frame()
lambda <- tmb_gh$lambda
gamma <- tmb_gh$gamma
delta <- tmb_gh$delta
if (BOOTSTRAP_SAMPLES != 0) {
  for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
    temp <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_TINN,
                                               mod = list(m = m,
                                                          lambda = tmb_gh$lambda,
                                                          gamma = tmb_gh$gamma,
                                                          delta = tmb_gh$delta),
                                               testing_params = list(m = m,
                                                                     lambda = tmb_gh$lambda,
                                                                     gamma = tmb_gh$gamma,
                                                                     delta = tmb_gh$delta))
    # The values from gamma are taken columnwise
    natural_parameters <- temp$natural_parameters
    natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
    bootstrap_tinn[idx_sample, 1:len_par] <- natural_parameters
  }
  if (BOOTSTRAP_SAMPLES == 1) {
    names(bootstrap_tinn) <- params_names_latex
  } else {
    colnames(bootstrap_tinn) <- params_names_latex
  }
  q <- apply(bootstrap_tinn,
             2,
             quantile.colwise)
  conf_int_tinn$Bootstrap.L <- q[1, ]
  conf_int_tinn$Bootstrap.U <- q[2, ]
}

# TMB confidence intervals --------------
# Manually cap values at their natural bound
# lambda must be strictly above 0
lambda_L <- pmax(0.0001,
                 tmb_gh$lambda - q95_norm * tmb_gh$lambda_std_error)
# gamma must be 0 or more
gamma_L <- pmax(0,
                tmb_gh$gamma - q95_norm * tmb_gh$gamma_std_error)
# delta must be 0 or above
delta_L <- pmax(0,
                tmb_gh$delta - q95_norm * tmb_gh$delta_std_error)
conf_int_tinn$TMB.L <- c(lambda_L,
                         gamma_L,
                         delta_L)
# no upper bound on lambda
# gamma must be 1 or less
gamma_U <- pmin(1,
                tmb_gh$gamma + q95_norm * tmb_gh$gamma_std_error)
# delta must be 1 or less
delta_U <- pmin(1,
                tmb_gh$delta + q95_norm * tmb_gh$delta_std_error)
conf_int_tinn$TMB.U <- c(tmb_gh$lambda + q95_norm * tmb_gh$lambda_std_error,
                         gamma_U,
                         delta_U)
# Coverage probabilities of the 3 CI methods -----------------
set.seed(FIRST_SEED + 4)
parameter_names <- paste0(rep("lambda",
                              m),
                          1:m)
for (col in 1:m) {
  # Get row and column indices for gamma instead of the default
  # columnwise index: the default indices are 1:m for the 1st column,
  # then (m + 1):(2 * m) for the 2nd, etc...
  parameter_names <- c(parameter_names,
                       paste0(sapply(X = 1:m,
                                     FUN = function(row) {paste0("gamma",
                                                                 row, ", ",
                                                                 col)})))
}
parameter_names <- c(parameter_names,
                     paste0(rep("delta",
                                m),
                            1:m))
# Record the times where the profile CI of a parameter successfully contains the parameter's true value
coverage_count_profile <-
  coverage_count_bootstrap <-
  coverage_count_tmb <-
  data.frame(parameter = parameter_names,
             count = 0,
             ratio = 0)
idx_coverage <- 0
while (idx_coverage < COVERAGE_SAMPLES) {
  idx_coverage <- idx_coverage + 1
  # Generate a data sample where nlminb converges
  # Loop as long as there is an issue with nlminb
  # Estimate a model
  coverage_model <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_TINN,
                                                       mod = list(m = m,
                                                                  lambda = tmb_gh$lambda,
                                                                  gamma = tmb_gh$gamma,
                                                                  delta = tmb_gh$delta),
                                                       testing_params = list(m = m,
                                                                             lambda = tmb_gh$lambda,
                                                                             gamma = tmb_gh$gamma,
                                                                             delta = tmb_gh$delta),
                                                       std_error = TRUE)
  
  # Save the occurrences of failures to generate a sample for which parameters can be estimated
  for (reason in c("state_number", "TMB_null", "TMB_converge", "TMB_G_null",
                   "TMB_G_converge", "TMB_H_null", "TMB_H_converge", "TMG_GH_null",
                   "TMG_GH_converge", "NA_value")) {
    coverage_skips_tinn[coverage_skips_tinn$m == m, reason] <-
      coverage_skips_tinn[coverage_skips_tinn$m == m, reason] +
      coverage_model$failure[reason]
  }
  
  # Confidence interval profiling -------------------------------------
  working_conf_int <- data.frame(w_parameter = w_params_names,
                                 lower = rep(NA,
                                             len_w_par),
                                 upper = rep(NA,
                                             len_w_par),
                                 stringsAsFactors = FALSE)
  for (idx_param in 1:len_w_par) {
    profile <- tmbprofile(obj = coverage_model$mod$obj,
                          name = idx_param,
                          trace = FALSE)
    
    ci <- tryCatch({
      confint(profile)
    },
    error = function(e){
      return(rep(NA,
                 2))
    })
    
    w_param_name <- w_params_names[idx_param]
    working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
    working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
  }
  
  lambda_profile_lower <- exp(working_conf_int$lower[lambda_indices])
  lambda_profile_upper <- exp(working_conf_int$upper[lambda_indices])
  gamma_profile_lower <- as.numeric(gamma.w2n(m,
                                              working_conf_int$lower[tgamma_indices]))
  gamma_profile_upper <- as.numeric(gamma.w2n(m,
                                              working_conf_int$upper[tgamma_indices]))
  
  # If profiling doesn't yield results for all parameters, try a new coverage sample
  estimates_coverage <- c(lambda_profile_lower,
                          lambda_profile_upper,
                          gamma_profile_lower,
                          gamma_profile_upper)
  test_null <- sapply(X = estimates_coverage, FUN = is.null)
  test_finite <- sapply(X = estimates_coverage, FUN = is.finite)
  # If some CI bounds are NULL or missing (NA) or infinite (Inf), try a new coverage sample
  if (any(test_null == TRUE) | any(test_finite == FALSE)) {
    idx_coverage <- idx_coverage - 1
    temp <- coverage_skips_tinn[coverage_skips_tinn$m == m, "profile"]
    coverage_skips_tinn[coverage_skips_tinn$m == m, "profile"] <- temp + 1
    next
  }
  # If the "true" value for lambda is in the CI, then increase the count
  real_lambda_profile_lower <- pmin(lambda_profile_lower,
                                    lambda_profile_upper)
  real_lambda_profile_upper <- pmax(lambda_profile_lower,
                                    lambda_profile_upper)
  indices <- which(tmb_gh$lambda >= real_lambda_profile_lower &
                     tmb_gh$lambda <= real_lambda_profile_upper)
  coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
  
  # Same for gamma
  real_gamma_profile_lower <- pmin(gamma_profile_lower,
                                   gamma_profile_upper)
  real_gamma_profile_upper <- pmax(gamma_profile_lower,
                                   gamma_profile_upper)
  indices <- which(as.vector(tmb_gh$gamma) >= real_gamma_profile_lower &
                     as.vector(tmb_gh$gamma) <= real_gamma_profile_upper)
  indices <- indices + m
  coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
  
  # Confidence interval bootstrap -----------------------------------
  bootstrap_tinn <- data.frame()
  if (BOOTSTRAP_SAMPLES != 0) {
    for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
      temp <- pois.HMM.generate.estimable.sample(
        ns = DATA_SIZE_TINN,
        mod = list(m = m,
                   lambda = coverage_model$natural_parameters$lambda,
                   gamma = coverage_model$natural_parameters$gamma,
                   delta = coverage_model$natural_parameters$delta),
        testing_params = list(m = m,
                              lambda = tmb_gh$lambda,
                              gamma = tmb_gh$gamma,
                              delta = tmb_gh$delta))
      # The values from gamma are taken columnwise
      natural_parameters <- temp$natural_parameters
      natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
      bootstrap_tinn[idx_sample, 1:len_par] <- natural_parameters
    }
    q <- apply(bootstrap_tinn,
               2,
               quantile.colwise)
    indices <- which(as.vector(tmb_gh$lambda) >= q[1, lambda_indices] &
                       as.vector(tmb_gh$lambda) <= q[2, lambda_indices])
    coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
    indices <- which(as.vector(tmb_gh$gamma) >= q[1, gamma_indices] &
                       as.vector(tmb_gh$gamma) <= q[2, gamma_indices]) + m
    coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
    indices <- which(as.vector(tmb_gh$delta) >= q[1, delta_indices] &
                       as.vector(tmb_gh$delta) <= q[2, delta_indices]) + m + m ^ 2
    coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
  }
  
  # Confidence interval TMB -----------------------------------------
  nat_par <- coverage_model$natural_parameters
  lambda_tmb_lower <- nat_par$lambda - q95_norm * nat_par$lambda_std_error
  lambda_tmb_upper <- nat_par$lambda + q95_norm * nat_par$lambda_std_error
  gamma_tmb_lower <- nat_par$gamma - q95_norm * nat_par$gamma_std_error
  gamma_tmb_upper <- nat_par$gamma + q95_norm * nat_par$gamma_std_error
  delta_tmb_lower <- nat_par$delta - q95_norm * nat_par$delta_std_error
  delta_tmb_upper <- nat_par$delta + q95_norm * nat_par$delta_std_error
  
  indices <- which(as.vector(tmb_gh$lambda) >= lambda_tmb_lower &
                     as.vector(tmb_gh$lambda) <= lambda_tmb_upper)
  coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
  indices <- which(as.vector(tmb_gh$gamma) >= gamma_tmb_lower &
                     as.vector(tmb_gh$gamma) <= gamma_tmb_upper)
  indices <- indices + m
  coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
  indices <- which(as.vector(tmb_gh$delta) >= delta_tmb_lower &
                     as.vector(tmb_gh$delta) <= delta_tmb_upper)
  indices <- indices + m + m ^ 2
  coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
  
}
coverage_count_profile[lambda_indices, "ratio"] <- coverage_count_profile[lambda_indices, "count"] /
  COVERAGE_SAMPLES
coverage_count_profile[gamma_indices, "ratio"] <- coverage_count_profile[gamma_indices, "count"] /
  COVERAGE_SAMPLES
# delta is not a parameter for us, so it has no profile CI
coverage_count_profile[delta_indices, "ratio"] <- NA

coverage_count_tmb$ratio <- coverage_count_tmb$count / COVERAGE_SAMPLES

coverage_count_bootstrap$ratio <- coverage_count_bootstrap$count / COVERAGE_SAMPLES

conf_int_tinn[1:(m ^ 2 + 2 * m), "Coverage.Profile"] <- coverage_count_profile$ratio * 100
conf_int_tinn[1:(m ^ 2 + 2 * m), "Coverage.Bootstrap"] <- coverage_count_bootstrap$ratio * 100
conf_int_tinn[1:(m ^ 2 + 2 * m), "Coverage.TMB"] <- coverage_count_tmb$ratio * 100

# Fixes -------------------------
# Fix label switching in conf_int_tinn
ordered_params <- pois.HMM.label.order(m = m,
                                       lambda = tmb_gh$lambda,
                                       gamma = tmb_gh$gamma,
                                       delta = tmb_gh$delta)

new_lambda_indices <- ordered_params$ordered_lambda_indices
new_gamma_indices <- ordered_params$ordered_gamma_vector_indices
new_delta_indices <- ordered_params$ordered_delta_indices

conf_int_tinn[lambda_indices, - 2] <- conf_int_tinn[lambda_indices, - 2][new_lambda_indices, ]
conf_int_tinn[gamma_indices, - 2] <- conf_int_tinn[gamma_indices, - 2][new_gamma_indices, ]
conf_int_tinn[delta_indices, - 2] <- conf_int_tinn[delta_indices, - 2][new_delta_indices, ]

# Reorder the TPM row-wise instead of column-wise
# Lexicographical parameter sort for gamma (sort on the parameter name)
new_gamma_indices_truncated_table <- order(conf_int_tinn[gamma_indices, "Parameter"])
# Replace rows by sorted rows
conf_int_tinn[gamma_indices, ] <- conf_int_tinn[gamma_indices, ][new_gamma_indices_truncated_table, ]

# The profile CIs may not be sorted, so we sort them manually
for (i in 1:length(conf_int_tinn[, 1])) {
  row <- conf_int_tinn[i, c("Profile.L", "Profile.U")]
  conf_int_tinn[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
}
conf_int_tinn$m <- as.integer(conf_int_tinn$m)

estim_benchmarks_df_tinn$m <- factor(estim_benchmarks_df_tinn$m,
                                     levels = M_LIST_TINN)
