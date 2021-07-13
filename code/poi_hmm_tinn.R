# Prepare the data & parameters, then estimate for different numbers of hidden states
for (idx in 1:length(M_LIST_TINN)) {
  set.seed(1)
  # Parameters and covariates --------------------------
  m <- M_LIST_TINN[idx]
  if (m == 1) {
    gamma_init <- matrix(1)
  } else {
    gamma_init <- matrix(0.2 / (m - 1), nrow = m, ncol = m)
    diag(gamma_init) <- 0.8
  }
  lambda_init <- seq(quantile(tinn_data, 0.1), quantile(tinn_data, 0.9), length.out = m)
  delta_init <- stat.dist(gamma_init)
  
  # Parameters & covariates for TMB ------------------
  working_params_init <- pois.HMM.pn2pw(m, lambda_init, gamma_init)
  TMB_data <- list(x = tinn_data, m = m)
  obj_init <- MakeADFun(TMB_data, working_params_init, DLL = "poi_hmm", silent = TRUE)
  parvect_init <- pois.HMM.pn2pw(m = m, lambda = lambda_init, gamma = gamma_init, delta = delta_init)
  parvect_init <- unlist(parvect_init)
  
  # Estimation ------------------------------------
  dm <- DM.estimate(x = tinn_data,
                    m = m,
                    lambda0 = lambda_init,
                    gamma0 = gamma_init)
  tmb <- TMB.estimate(TMB_data = TMB_data,
                      parameters = working_params_init,
                      MakeADFun_obj = obj_init)
  tmb_g <- TMB.estimate(TMB_data = TMB_data,
                        parameters = working_params_init,
                        MakeADFun_obj = obj_init,
                        gradient = TRUE)
  tmb_h <- TMB.estimate(TMB_data = TMB_data,
                        parameters = working_params_init,
                        MakeADFun_obj = obj_init,
                        hessian = TRUE)
  tmb_gh <- TMB.estimate(TMB_data = TMB_data,
                         parameters = working_params_init,
                         MakeADFun_obj = obj_init,
                         gradient = TRUE,
                         hessian = TRUE,
                         std_error = TRUE)
  
  # If one doesn't converge successfully, stop
  if (dm$convergence != 0) {
    stop(paste("dm didn't converge properly, tinn dataset, m =", m))
  }
  if (tmb$convergence != 0) {
    stop(paste("tmb didn't converge properly, tinn dataset, m =", m))
  }
  if (tmb_g$convergence != 0) {
    stop(paste("tmb_g didn't converge properly, tinn dataset, m =", m))
  }
  if (tmb_h$convergence != 0) {
    stop(paste("tmb_h didn't converge properly, tinn dataset, m =", m))
  }
  if (tmb_gh$convergence != 0) {
    stop(paste("tmb_gh didn't converge properly, tinn dataset, m =", m))
  }
  
  # Store negative log-likelihoods --------------------------------------------
  # We want to store the next 5 values without using rbind
  # If the file is rerun, this variable won't store new results on top of the
  # old ones
  indices <- (1 + 5 * (idx - 1)):(5 + 5 * (idx - 1))
  mllk_values_tinn[indices, ] <- cbind(m = rep(m, 5),
                                       procedure = PROCEDURES,
                                       mllk = c(dm$mllk,
                                                tmb$mllk,
                                                tmb_g$mllk,
                                                tmb_h$mllk,
                                                tmb_gh$mllk),
                                       AIC = c(dm$AIC,
                                               tmb$AIC,
                                               tmb_g$AIC,
                                               tmb_h$AIC,
                                               tmb_gh$AIC),
                                       BIC = c(dm$BIC,
                                               tmb$BIC,
                                               tmb_g$BIC,
                                               tmb_h$BIC,
                                               tmb_gh$BIC))
  
  # Creating variables for the CIs -----------------
  tmb_CI <- tmb_gh
  params_names_latex <- paste0(rep("$\\lambda_{", m), 1:m, "}$")
  for (gamma_idx in 1:m ^ 2) {
    # Get row and column indices for gamma instead of the default
    # columnwise index: the default indices are 1:m for the 1st column,
    # then (m + 1):(2 * m) for the 2nd, etc...
    row_col_idx <- matrix.col.idx.to.rowcol(gamma_idx, m)
    params_names_latex <- c(params_names_latex,
                            paste0("$\\gamma_{", paste0(row_col_idx, collapse = ""), "}$"))
  }
  params_names_latex <- c(params_names_latex,
                          paste0(rep("$\\delta_{", m), 1:m, "}$"))
  len_par <- length(params_names_latex)
  indices <- (length(conf_int_tinn$m) + 1):(length(conf_int_tinn$m) + len_par)
  conf_int_tinn[indices, "m"] <- m
  conf_int_tinn[indices, "Parameter"] <- params_names_latex
  # Reminder, PARAMS_NAMES contains c("lambda", "gamma", "delta")
  conf_int_tinn[indices, "Estimate"] <- unlist(tmb_CI[PARAMS_NAMES])
  
  param_tmb_CI <- pois.HMM.pn2pw(m = m, lambda = tmb_CI$lambda, gamma = tmb_CI$gamma)
  
  if (m == 1) {
    w_params_names <- c("tlambda1")
  } else {
    w_params_names <- c(paste0(rep("tlambda", m), 1:m),
                        paste0(rep("tgamma", m ^ 2 - m), 1:(m ^ 2 - m)))
  }
  len_w_par <- length(w_params_names)
  working_conf_int <- data.frame(w_parameter = w_params_names,
                                 lower = rep(NA, len_w_par),
                                 upper = rep(NA, len_w_par),
                                 stringsAsFactors = FALSE)
  
  lambda_indices <- 1:m
  gamma_indices <- m + 1:(m ^ 2)
  delta_indices <- m ^ 2 + m + (1:m)
  tgamma_indices <- (m + 1):(m ^ 2)
  
  # Benchmark the same dataset many times to check the benchmark durations has low variance --------------
  set.seed(2)
  if (CONSISTENCY_BENCHMARK_TINN != 0) {
    # Parameters & covariates for DM and TMB are already done
    # Estimation benchmark
    temp <- microbenchmark("DM" = nlminb(parvect_init, pois.HMM.mllk, x_alias = tinn_data, m_alias = m)$convergence==0,
                           "TMB" = nlminb(obj_init$par, obj_init$fn)$convergence==0,
                           "TMB_G" = nlminb(obj_init$par, obj_init$fn, gradient = obj_init$gr)$convergence==0,
                           "TMB_H" = nlminb(obj_init$par, obj_init$fn, hessian = obj_init$he)$convergence==0,
                           "TMB_GH" = nlminb(obj_init$par, obj_init$fn, gradient = obj_init$gr, hessian = obj_init$he)$convergence==0,
                           times = CONSISTENCY_BENCHMARK_TINN,
                           setup = obj_init <<- MakeADFun(TMB_data, working_params_init, DLL = "poi_hmm", silent = TRUE),
                           check = "equal")
    
    iterations <- c(nlminb(parvect_init, pois.HMM.mllk, x_alias = tinn_data, m_alias = m)$iterations,
                    nlminb(obj_init$par, obj_init$fn)$iteration,
                    nlminb(obj_init$par, obj_init$fn, gradient = obj_init$gr)$iteration,
                    nlminb(obj_init$par, obj_init$fn, hessian = obj_init$he)$iteration,
                    nlminb(obj_init$par, obj_init$fn, gradient = obj_init$gr, hessian = obj_init$he)$iteration)
    
    times <- temp$time/10^9
    timeDM <- times[temp$expr == "DM"]
    timeTMB <- times[temp$expr == "TMB"]
    timeTMB_G <- times[temp$expr == "TMB_G"]
    timeTMB_H <- times[temp$expr == "TMB_H"]
    timeTMB_GH <- times[temp$expr == "TMB_GH"]
    consistency_estim_benchmarks_df_tinn <- rbind(consistency_estim_benchmarks_df_tinn,
                                                  data.frame(time = c(timeDM, timeTMB, timeTMB_G, timeTMB_H, timeTMB_GH),
                                                             m = rep(m, length(PROCEDURES) * CONSISTENCY_BENCHMARK_TINN),
                                                             procedure = rep(PROCEDURES, each = CONSISTENCY_BENCHMARK_TINN),
                                                             iterations = rep(iterations, each = CONSISTENCY_BENCHMARK_TINN)))
  }
  
  # Benchmarks ------------
  set.seed(3)
  if (BENCHMARK_SAMPLES != 0) {
    for (idx_counter in 1:BENCHMARK_SAMPLES) {
      # Generate data that can be estimated by TMB_GH
      # and is tested on the slightly off parameters from the beginning of this file
      # The goal is to have a dataset that poses no estimation problem,
      # when estimated with guessed initial parameters
      benchmark_model <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_TINN,
                                                            mod = list(m = m,
                                                                       lambda = tmb_CI$lambda,
                                                                       gamma = tmb_CI$gamma,
                                                                       delta = tmb_CI$delta),
                                                            testing_params = list(m = m,
                                                                                  lambda = lambda_init,
                                                                                  gamma = gamma_init,
                                                                                  delta = delta_init),
                                                            test_marqLevAlg = TRUE)
      benchmark_data <- benchmark_model$data
      # Benchmark all different combinations of gradient and hessians with DM ----------------
      # Parameters & covariates for DM and TMB
      TMB_benchmark_data <- list(x = benchmark_data, m = m)
      obj_benchmark <- MakeADFun(TMB_benchmark_data, working_params_init, DLL = "poi_hmm", silent = TRUE)
      parvect_benchmark <- pois.HMM.pn2pw(m = m, lambda = lambda_init, gamma = gamma_init, delta = delta_init)
      # nlminb needs a vector, not a list
      parvect_benchmark <- unlist(parvect_benchmark)
      # Estimation benchmark
      temp <- microbenchmark("DM" = nlminb(parvect_benchmark, pois.HMM.mllk, x_alias = benchmark_data, m_alias = m)$convergence==0,
                             "TMB" = nlminb(obj_benchmark$par, obj_benchmark$fn)$convergence==0,
                             "TMB_G" = nlminb(obj_benchmark$par, obj_benchmark$fn, gradient = obj_benchmark$gr)$convergence==0,
                             "TMB_H" = nlminb(obj_benchmark$par, obj_benchmark$fn, hessian = obj_benchmark$he)$convergence==0,
                             "TMB_GH" = nlminb(obj_benchmark$par, obj_benchmark$fn, gradient = obj_benchmark$gr, hessian = obj_benchmark$he)$convergence==0,
                             times = 1,
                             check = "equal",
                             setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data, working_params_init, DLL = "poi_hmm", silent = TRUE))
      times <- temp$time/10^9
      timeDM <- times[temp$expr == "DM"]
      timeTMB <- times[temp$expr == "TMB"]
      timeTMB_G <- times[temp$expr == "TMB_G"]
      timeTMB_H <- times[temp$expr == "TMB_H"]
      timeTMB_GH <- times[temp$expr == "TMB_GH"]
      
      iterDM <- nlminb(parvect_benchmark, pois.HMM.mllk, x_alias = benchmark_data, m_alias = m)$iterations
      iterTMB <- nlminb(obj_benchmark$par, obj_benchmark$fn)$iterations
      iterTMB_G <- nlminb(obj_benchmark$par, obj_benchmark$fn, gradient = obj_benchmark$gr)$iterations
      iterTMB_H <- nlminb(obj_benchmark$par, obj_benchmark$fn, hessian = obj_benchmark$he)$iterations
      iterTMB_GH <- nlminb(obj_benchmark$par, obj_benchmark$fn, gradient = obj_benchmark$gr, hessian = obj_benchmark$he)$iterations
      
      estim_benchmarks_df_tinn <- rbind(estim_benchmarks_df_tinn,
                                        data.frame(time = c(timeDM, timeTMB, timeTMB_G, timeTMB_H, timeTMB_GH),
                                                   m = rep(m, length(PROCEDURES)),
                                                   procedure = PROCEDURES,
                                                   iterations = c(iterDM, iterTMB, iterTMB_G, iterTMB_H, iterTMB_GH),
                                                   dataset_number = rep(idx_counter, length(PROCEDURES))))
      
      # Benchmark mllk times --------------------------------------------------
      # tmb_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
      #                               parameters = parvect_benchmark,
      #                               MakeADFun_obj = obj_benchmark)
      # tmb_g_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
      #                                 parameters = parvect_benchmark,
      #                                 MakeADFun_obj = obj_benchmark,
      #                                 gradient = TRUE)
      # tmb_h_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
      #                                 parameters = parvect_benchmark,
      #                                 MakeADFun_obj = obj_benchmark,
      #                                 hessian = TRUE)
      tmb_gh_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
                                       parameters = parvect_benchmark,
                                       MakeADFun_obj = obj_benchmark,
                                       gradient = TRUE,
                                       hessian = TRUE)
      
      # param_tmb <- pois.HMM.pn2pw(m = m, lambda = tmb_benchmark$lambda, gamma = tmb_benchmark$gamma)
      # param_tmb_g <- pois.HMM.pn2pw(m = m, lambda = tmb_g_benchmark$lambda, gamma = tmb_g_benchmark$gamma)
      # param_tmb_h <- pois.HMM.pn2pw(m = m, lambda = tmb_h_benchmark$lambda, gamma = tmb_h_benchmark$gamma)
      param_tmb_gh <- pois.HMM.pn2pw(m = m, lambda = tmb_gh_benchmark$lambda, gamma = tmb_gh_benchmark$gamma)
      # model1 <- MakeADFun(TMB_data, param_tmb, DLL = "poi_hmm", silent = TRUE)
      # model2 <- MakeADFun(TMB_data, param_tmb_g, DLL = "poi_hmm", silent = TRUE)
      # model3 <- MakeADFun(TMB_data, param_tmb_h, DLL = "poi_hmm", silent = TRUE)
      model4 <- MakeADFun(TMB_data, param_tmb_gh, DLL = "poi_hmm", silent = TRUE)
      parvect_mllk <- pois.HMM.pn2pw(m, lambda = tmb_gh_benchmark$lambda, gamma = tmb_gh_benchmark$gamma)
      
      temp <- microbenchmark("DM" = pois.HMM.mllk(parvect_mllk, tinn_data, m),
                             # "TMB" = model1$fn(model1$par),
                             # "TMB_G" = model2$fn(model2$par),
                             # "TMB_H" = model3$fn(model3$par),
                             "TMB_GH" = model4$fn(model4$par),
                             times = 1)
      
      times <- temp$time/10^9
      timeDM <- times[temp$expr == "DM"]
      # timeTMB <- times[temp$expr == "TMB"]
      # timeTMB_G <- times[temp$expr == "TMB_G"]
      # timeTMB_H <- times[temp$expr == "TMB_H"]
      timeTMB_GH <- times[temp$expr == "TMB_GH"]
      mllk_times_df_tinn <- rbind(mllk_times_df_tinn,
                                  data.frame(time = c(timeDM, timeTMB_GH),
                                             m = rep(m, 2),
                                             procedure = PROCEDURES[c(1,5)],
                                             dataset_number = rep(idx_counter, 2)))
      
      # Benchmark different optimization methods ----------------------------------------------
      temp <- microbenchmark("BFGS" = optim(par = obj_benchmark$par, fn = obj_benchmark$fn, gr = obj_benchmark$gr, method = "BFGS", control = ctrl)$convergence==0,
                             "L-BFGS-B" = optim(par = obj_benchmark$par, fn = obj_benchmark$fn, gr = obj_benchmark$gr, method = "L-BFGS-B", control = ctrl)$convergence==0,
                             "nlm" = nlm(f = nlmfn, p = obj_benchmark$par, obj_benchmark, iterlim = 10000)$code==1,
                             "nlminb" = nlminb(start = obj_benchmark$par, objective = obj_benchmark$fn, gradient = obj_benchmark$gr, hessian = obj_benchmark$he)$convergence==0,
                             "hjn" = hjn(par = obj_benchmark$par, fn = obj_benchmark$fn)$convergence==0,
                             "marqLevAlg" = marqLevAlg(b = obj_benchmark$par, fn = obj_benchmark$fn, gr = obj_benchmark$gr, hess = obj_benchmark$he, maxiter = 10000)$istop==1,
                             times = 1,
                             check = "equal",
                             setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data, working_params_init, DLL = "poi_hmm", silent = TRUE))
      
      times <- temp$time/10^9
      timeBFGS <- times[temp$expr == "BFGS"]
      timeL_BFGS_B <- times[temp$expr == "L-BFGS-B"]
      timenlm <- times[temp$expr == "nlm"]
      timenlminb <- times[temp$expr == "nlminb"]
      timehjn <- times[temp$expr == "hjn"]
      timemarqLevAlg <- times[temp$expr == "marqLevAlg"]
      method_comparison_df_tinn <- rbind(method_comparison_df_tinn,
                                         data.frame(time = c(timeBFGS, timeL_BFGS_B, timenlm, timenlminb, timehjn, timemarqLevAlg),
                                                    m = rep(m, length(PROCEDURES_METHOD)),
                                                    procedure = PROCEDURES_METHOD,
                                                    dataset_number = idx_counter))
    }
  }
  
  # Profiling the likelihood --------------------------
  for (idx_param in 1:len_w_par) {
    profile <- tmbprofile(obj = tmb_CI$obj,
                          name = idx_param,
                          trace = FALSE)
    
    ci <- tryCatch({
      confint(profile)
    },
    error = function(e){
      return(rep(NA, 2))
    })
    
    w_param_name <- w_params_names[idx_param]
    working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
    working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
  }
  
  # Transform the working parameters into natural ones
  # Lambda (m values)
  conf_int_tinn$Profile.L[which(conf_int_tinn$m == m)][lambda_indices] <- exp(working_conf_int$lower[lambda_indices])
  # Gamma (m^2-m working parameters, m^2 natural ones)
  if (!anyNA(working_conf_int$lower[tgamma_indices])) {
    natural_gamma <- as.numeric(gamma.w2n(m, working_conf_int$lower[tgamma_indices]))
    conf_int_tinn$Profile.L[which(conf_int_tinn$m == m)][gamma_indices] <- natural_gamma
  }
  # Lambda (m values)
  conf_int_tinn$Profile.U[which(conf_int_tinn$m == m)][lambda_indices] <- exp(working_conf_int$upper[lambda_indices])
  # Gamma (m^2-m working parameters, m^2 natural ones)
  if (!anyNA(working_conf_int$upper[tgamma_indices])) {
    natural_gamma <- as.numeric(gamma.w2n(m, working_conf_int$upper[tgamma_indices]))
    conf_int_tinn$Profile.U[which(conf_int_tinn$m == m)][gamma_indices] <- natural_gamma
  }
  
  # Bootstrap ---------------------------
  set.seed(4)
  bootstrap_tinn <- data.frame()
  lambda <- tmb_CI$lambda
  gamma <- tmb_CI$gamma
  delta <- tmb_CI$delta
  if (BOOTSTRAP_SAMPLES != 0) {
    for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
      temp <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_TINN,
                                                 mod = list(m = m,
                                                            lambda = tmb_CI$lambda,
                                                            gamma = tmb_CI$gamma,
                                                            delta = tmb_CI$delta),
                                                 testing_params = list(m = m,
                                                                       lambda = lambda_init,
                                                                       gamma = gamma_init,
                                                                       delta = delta_init))$natural_parameters
      # The values from gamma are taken columnwise
      natural_parameters <- temp
      natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
      bootstrap_tinn[idx_sample, 1:len_par] <- natural_parameters
    }
    colnames(bootstrap_tinn) <- params_names_latex
    q <- apply(bootstrap_tinn, 2, quantile.colwise)
    conf_int_tinn$Bootstrap.L[which(conf_int_tinn$m == m)] <- q[1, ]
    conf_int_tinn$Bootstrap.U[which(conf_int_tinn$m == m)] <- q[2, ]
  }
  
  # TMB confidence intervals --------------
  # Manually cap values at their natural bound
  # lambda must be strictly above 0
  lambda_L <- pmax(0.0001, tmb_gh$lambda - q95_norm * tmb_gh$lambda_std_error)
  # gamma must be 0 or more
  gamma_L <- pmax(0, tmb_gh$gamma - q95_norm * tmb_gh$gamma_std_error)
  # delta must be 0 or above
  delta_L <- pmax(0, tmb_gh$delta - q95_norm * tmb_gh$delta_std_error)
  conf_int_tinn$TMB.L[which(conf_int_tinn$m == m)] <- c(lambda_L,
                                                        gamma_L,
                                                        delta_L)
  # no upper bound on lambda
  # gamma must be 1 or less
  gamma_U <- pmin(1, tmb_gh$gamma + q95_norm * tmb_gh$gamma_std_error)
  # delta must be 1 or less
  delta_U <- pmin(1, tmb_gh$delta + q95_norm * tmb_gh$delta_std_error)
  conf_int_tinn$TMB.U[which(conf_int_tinn$m == m)] <- c(tmb_gh$lambda + q95_norm * tmb_gh$lambda_std_error,
                                                        gamma_U,
                                                        delta_U)
  # Coverage probabilities of the 3 CI methods -----------------
  set.seed(5)
  parameter_names <- paste0(rep("lambda", m), 1:m)
  for (gamma_idx in 1:(m ^ 2)) {
    # Get row and column indices for gamma instead of the default
    # columnwise index: the default indices are 1:m for the 1st column,
    # then (m + 1):(2 * m) for the 2nd, etc...
    row_col_idx <- matrix.col.idx.to.rowcol(gamma_idx, m)
    parameter_names <- c(parameter_names,
                         paste0("gamma", toString(row_col_idx)))
  }
  parameter_names <- c(parameter_names, paste0(rep("delta", m), 1:m))
  coverage_count_profile <- coverage_count_bootstrap <- coverage_count_tmb <- data.frame(parameter = parameter_names,
                                                                                         count = 0,
                                                                                         ratio = 0)
  idx_coverage <- 0
  while (idx_coverage < COVERAGE_SAMPLES) {
    idx_coverage <- idx_coverage + 1
    # Generate a data sample where nlminb converges
    # Loop as long as there is an issue with nlminb
    # Estimate a model
    coverage_model <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_TINN,
                                                         mod = list(m = m,
                                                                    lambda = tmb_CI$lambda,
                                                                    gamma = tmb_CI$gamma,
                                                                    delta = tmb_CI$delta),
                                                         testing_params = list(m = m,
                                                                               lambda = lambda_init,
                                                                               gamma = gamma_init,
                                                                               delta = delta_init),
                                                         std_error = TRUE)
    
    # Save the occurrences of failures to generate a sample for which parameters can be estimated
    for (reason in c("state_number", "TMB_null", "TMB_converge", "TMB_G_null",
                     "TMB_G_converge", "TMB_H_null", "TMB_H_converge", "TMG_GH_null",
                     "TMG_GH_converge", "marqLevAlg_converge", "NA_value")) {
      coverage_skips_tinn[coverage_skips_tinn$m == m, reason] <- coverage_skips_tinn[coverage_skips_tinn$m == m, reason] + coverage_model$failure[reason]
    }
    
    # Confidence interval profiling -------------------------------------
    working_conf_int <- data.frame(w_parameter = w_params_names,
                                   lower = rep(NA, len_w_par),
                                   upper = rep(NA, len_w_par),
                                   stringsAsFactors = FALSE)
    for (idx_param in 1:len_w_par) {
      profile <- tmbprofile(obj = coverage_model$mod$obj,
                            name = idx_param,
                            trace = FALSE)
      
      ci <- tryCatch({
        confint(profile)
      },
      error = function(e){
        return(rep(NA, 2))
      })
      
      w_param_name <- w_params_names[idx_param]
      working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
      working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
    }
    
    lambda_profile_lower <- exp(working_conf_int$lower[lambda_indices])
    lambda_profile_upper <- exp(working_conf_int$upper[lambda_indices])
    gamma_profile_lower <- as.numeric(gamma.w2n(m, working_conf_int$lower[tgamma_indices]))
    gamma_profile_upper <- as.numeric(gamma.w2n(m, working_conf_int$upper[tgamma_indices]))
    
    # If profiling doesn't yield results for all parameters, try a new coverage sample
    if (anyNA(c(lambda_profile_lower, lambda_profile_upper, gamma_profile_lower, gamma_profile_upper), recursive = TRUE)) {
      idx_coverage <- idx_coverage - 1
      coverage_skips_tinn[coverage_skips_tinn$m == m, "profile"] <- coverage_skips_tinn[coverage_skips_tinn$m == m, "profile"] + 1
      next
    }
    # If the "true" value for lambda is in the CI, then increase the count
    real_lambda_profile_lower <- pmin(lambda_profile_lower, lambda_profile_upper)
    real_lambda_profile_upper <- pmax(lambda_profile_lower, lambda_profile_upper)
    indices <- which(tmb_CI$lambda >= real_lambda_profile_lower & tmb_CI$lambda <= real_lambda_profile_upper)
    coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
    
    # Same for gamma
    real_gamma_profile_lower <- pmin(gamma_profile_lower, gamma_profile_upper)
    real_gamma_profile_upper <- pmax(gamma_profile_lower, gamma_profile_upper)
    indices <- which(as.vector(tmb_CI$gamma) >= real_gamma_profile_lower & as.vector(tmb_CI$gamma) <= real_gamma_profile_upper)
    indices <- indices + m
    coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
    
    # # Sometimes, there are no profile CIs for gamma because of the need to transform the CIs to their natural form
    # # In these cases, we reduce the total coverage_adjusted_sample_size_gamma_tinn[idx, "size"]
    # # Either or both bounds for gamma can be missing due to missing tgamma values
    # # If there are missing tgamma values, the confidence interval has 0 or 1 bound
    # # and is thus useless for coverage probabilities
    # if (anyNA(gamma_profile_lower, recursive = TRUE) || anyNA(gamma_profile_upper, recursive = TRUE)) {
    #   coverage_adjusted_sample_size_gamma_tinn[idx, "size"] <- coverage_adjusted_sample_size_gamma_tinn[idx, "size"] - 1
    # } else {
    #   # Same as above for gamma
    #   indices <- which(as.vector(tmb_CI$gamma) >= gamma_profile_lower & as.vector(tmb_CI$gamma) <= gamma_profile_upper)
    #   indices <- indices + m
    #   coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
    # }
    
    # Confidence interval bootstrap -----------------------------------
    bootstrap_tinn <- data.frame()
    if (BOOTSTRAP_SAMPLES != 0) {
      for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
        temp <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_TINN,
                                                   mod = list(m = m,
                                                              lambda = coverage_model$mod$lambda,
                                                              gamma = coverage_model$mod$gamma,
                                                              delta = coverage_model$mod$delta),
                                                   testing_params = list(m = m,
                                                                         lambda = tmb_CI$lambda,
                                                                         gamma = tmb_CI$gamma,
                                                                         delta = tmb_CI$delta))$natural_parameters
        # The values from gamma are taken columnwise
        natural_parameters <- temp
        natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
        bootstrap_tinn[idx_sample, 1:len_par] <- natural_parameters
      }
      # colnames(bootstrap_tinn) <- params_names_latex
      q <- apply(bootstrap_tinn, 2, quantile.colwise)
      indices <- which(as.vector(tmb_CI$lambda) >= q[1, lambda_indices] & as.vector(tmb_CI$lambda) <= q[2, lambda_indices])
      coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
      indices <- which(as.vector(tmb_CI$gamma) >= q[1, gamma_indices] & as.vector(tmb_CI$gamma) <= q[2, gamma_indices]) + m
      coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
      indices <- which(as.vector(tmb_CI$delta) >= q[1, delta_indices] & as.vector(tmb_CI$delta) <= q[2, delta_indices]) + m + m ^ 2
      coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
    }
    
    # Confidence interval TMB -----------------------------------------
    lambda_tmb_lower <- coverage_model$natural_parameters$lambda - q95_norm * coverage_model$natural_parameters$lambda_std_error
    lambda_tmb_upper <- coverage_model$natural_parameters$lambda + q95_norm * coverage_model$natural_parameters$lambda_std_error
    gamma_tmb_lower <- coverage_model$natural_parameters$gamma - q95_norm * coverage_model$natural_parameters$gamma_std_error
    gamma_tmb_upper <- coverage_model$natural_parameters$gamma + q95_norm * coverage_model$natural_parameters$gamma_std_error
    delta_tmb_lower <- coverage_model$natural_parameters$delta - q95_norm * coverage_model$natural_parameters$delta_std_error
    delta_tmb_upper <- coverage_model$natural_parameters$delta + q95_norm * coverage_model$natural_parameters$delta_std_error
    
    indices <- which(as.vector(tmb_CI$lambda) >= lambda_tmb_lower & as.vector(tmb_CI$lambda) <= lambda_tmb_upper)
    coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
    indices <- which(as.vector(tmb_CI$gamma) >= gamma_tmb_lower & as.vector(tmb_CI$gamma) <= gamma_tmb_upper)
    indices <- indices + m
    coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
    indices <- which(as.vector(tmb_CI$delta) >= delta_tmb_lower & as.vector(tmb_CI$delta) <= delta_tmb_upper)
    indices <- indices + m + m ^ 2
    coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
    
  }
  coverage_count_profile[lambda_indices, "ratio"] <- coverage_count_profile[lambda_indices, "count"] / COVERAGE_SAMPLES
  # coverage_count_profile[gamma_indices, "ratio"] <- coverage_count_profile[gamma_indices, "count"] / coverage_adjusted_sample_size_gamma_tinn[idx, "size"]
  coverage_count_profile[gamma_indices, "ratio"] <- coverage_count_profile[gamma_indices, "count"] / COVERAGE_SAMPLES
  coverage_count_profile[delta_indices, "ratio"] <- NA # delta is not a parameter for us, so it has no profile CI
  
  coverage_count_tmb$ratio <- coverage_count_tmb$count / COVERAGE_SAMPLES
  
  coverage_count_bootstrap$ratio <- coverage_count_bootstrap$count / COVERAGE_SAMPLES
  
  conf_int_tinn[conf_int_tinn$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Profile"] <- coverage_count_profile$ratio * 100
  conf_int_tinn[conf_int_tinn$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Bootstrap"] <- coverage_count_bootstrap$ratio * 100
  conf_int_tinn[conf_int_tinn$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.TMB"] <- coverage_count_tmb$ratio * 100
  
  # Fixes -------------------------
  # Fix label switching in conf_int_tinn
  ordered_params <- pois.HMM.label.order(m = m,
                                         lambda = tmb_CI$lambda,
                                         gamma = tmb_CI$gamma,
                                         delta = tmb_CI$delta)
  
  new_lambda_indices <- ordered_params$ordered_lambda_indices
  new_gamma_indices <- ordered_params$ordered_gamma_vector_indices
  new_delta_indices <- ordered_params$ordered_delta_indices
  
  conf_int_tinn[conf_int_tinn$m == m, - 2][lambda_indices, ] <- conf_int_tinn[conf_int_tinn$m == m, - 2][lambda_indices, ][new_lambda_indices, ]
  conf_int_tinn[conf_int_tinn$m == m, - 2][gamma_indices, ] <- conf_int_tinn[conf_int_tinn$m == m, - 2][gamma_indices, ][new_gamma_indices, ]
  conf_int_tinn[conf_int_tinn$m == m, - 2][delta_indices, ] <- conf_int_tinn[conf_int_tinn$m == m, - 2][delta_indices, ][new_delta_indices, ]
  
  # Reorder the TPM row-wise instead of column-wise
  # Lexicographical parameter sort for gamma (sort on the parameter name)
  new_gamma_indices_truncated_table <- order(conf_int_tinn[conf_int_tinn$m == m, ][gamma_indices, "Parameter"])
  # Replace rows by sorted rows
  conf_int_tinn[conf_int_tinn$m == m, ][gamma_indices, ] <- conf_int_tinn[conf_int_tinn$m == m, ][gamma_indices, ][new_gamma_indices_truncated_table, ]
}

# The profile CIs may not be sorted, so we sort them manually
for (i in 1:length(conf_int_tinn[, 1])) {
  row <- conf_int_tinn[i, c("Profile.L", "Profile.U")]
  conf_int_tinn[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
}
mllk_values_tinn$m <- as.factor(as.numeric(mllk_values_tinn$m))
mllk_values_tinn$mllk <- as.numeric(mllk_values_tinn$mllk)
mllk_values_tinn$AIC <- as.numeric(mllk_values_tinn$AIC)
mllk_values_tinn$BIC <- as.numeric(mllk_values_tinn$BIC)
conf_int_tinn$m <- as.integer(conf_int_tinn$m)

# Remove gamma and delta when m = 1 because they're pointless
if (1 %in% M_LIST_TINN) {
  conf_int_tinn <- conf_int_tinn[- c(2,3), ]
}

estim_benchmarks_df_tinn$m <- factor(estim_benchmarks_df_tinn$m, levels = M_LIST_TINN)