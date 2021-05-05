# Prepare the data & parameters, then estimate for different numbers of hidden states
for (idx in 1:length(M_LIST_LAMB)) {
  # Parameters and covariates --------------------------
  m <- M_LIST_LAMB[idx]
  if (m == 1) {
    gamma_init <- matrix(1)
  } else {
    gamma_init <- matrix(0.2 / (m - 1), nrow = m, ncol = m)
    diag(gamma_init) <- 0.8
  }
  lambda_init <- seq(0.3, 4, length.out = m)
  delta_init <- stat.dist(gamma_init)
  
  # Parameters & covariates for TMB ------------------
  working_params_init <- pois.HMM.pn2pw(m, lambda_init, gamma_init)
  TMB_data <- list(x = lamb_data, m = m)
  obj_init <- MakeADFun(TMB_data, working_params_init, DLL = "poi_hmm", silent = TRUE)
  parvect_init <- pois.HMM.pn2pw(m = m, lambda = lambda_init, gamma = gamma_init, delta = delta_init)
  parvect_init <- unlist(parvect_init)
  
  # Estimation ------------------------------------
  dm <- DM.estimate(x = lamb_data,
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
    stop(paste("dm didn't converge properly, lamb dataset, m =", m))
  }
  if (tmb$convergence != 0) {
    stop(paste("tmb didn't converge properly, lamb dataset, m =", m))
  }
  if (tmb_g$convergence != 0) {
    stop(paste("tmb_g didn't converge properly, lamb dataset, m =", m))
  }
  if (tmb_h$convergence != 0) {
    stop(paste("tmb_h didn't converge properly, lamb dataset, m =", m))
  }
  if (tmb_gh$convergence != 0) {
    stop(paste("tmb_gh didn't converge properly, lamb dataset, m =", m))
  }
  
  # Store negative log-likelihoods --------------------------------------------
  # We want to store the next 5 values without using rbind
  # If the file is rerun, this variable won't store new results on top of the
  # old ones
  indexes <- (1 + 5 * (idx - 1)):(5 + 5 * (idx - 1))
  mllk_values_lamb[indexes, ] <- cbind(m = rep(m, 5),
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
    # Get row and column indexes for gamma instead of the default
    # columnwise index: the default indexes are 1:m for the 1st column,
    # then (m + 1):(2 * m) for the 2nd, etc...
    row_col_idx <- matrix.col.idx.to.rowcol(gamma_idx, m)
    params_names_latex <- c(params_names_latex,
                            paste0("$\\gamma_{", toString(row_col_idx), "}$"))
  }
  params_names_latex <- c(params_names_latex,
                          paste0(rep("$\\delta_{", m), 1:m, "}$"))
  len_par <- length(params_names_latex)
  indexes <- (length(conf_int_lamb$m) + 1):(length(conf_int_lamb$m) + len_par)
  conf_int_lamb[indexes, "m"] <- m
  conf_int_lamb[indexes, "Parameter"] <- params_names_latex
  # Reminder, PARAMS_NAMES contains c("lambda", "gamma", "delta")
  conf_int_lamb[indexes, "Estimate"] <- unlist(tmb_CI[PARAMS_NAMES])
  
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
  
  # Benchmark the same dataset many times to check the benchmark durations has low variance --------------
  # Parameters & covariates for DM and TMB are already done
  # Estimation benchmark
  temp <- microbenchmark("DM" = nlminb(parvect_init, pois.HMM.mllk, x_alias = lamb_data, m_alias = m)$convergence==0,
                         "TMB" = nlminb(obj_init$par, obj_init$fn)$convergence==0,
                         "TMB_G" = nlminb(obj_init$par, obj_init$fn, gradient = obj_init$gr)$convergence==0,
                         "TMB_H" = nlminb(obj_init$par, obj_init$fn, hessian = obj_init$he)$convergence==0,
                         "TMB_GH" = nlminb(obj_init$par, obj_init$fn, gradient = obj_init$gr, hessian = obj_init$he)$convergence==0,
                         times = CONSISTENCY_BENCHMARK_LAMB,
                         control = list(order = "block", warmup = CONSISTENCY_BENCHMARK_LAMB / 4),
                         setup = obj_init <<- MakeADFun(TMB_data, working_params_init, DLL = "poi_hmm", silent = TRUE),
                         check = "equal")
  iterations <- c(nlminb(parvect_init, pois.HMM.mllk, x_alias = lamb_data, m_alias = m)$iterations,
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
  consistency_estim_benchmarks_df_lamb <- rbind(consistency_estim_benchmarks_df_lamb,
                                                data.frame(time = c(timeDM, timeTMB, timeTMB_G, timeTMB_H, timeTMB_GH),
                                                           m = rep(m, length(PROCEDURES) * CONSISTENCY_BENCHMARK_LAMB),
                                                           procedure = rep(PROCEDURES, each = CONSISTENCY_BENCHMARK_LAMB),
                                                           iterations = rep(iterations, each = CONSISTENCY_BENCHMARK_LAMB)))
  
  # Benchmarks ------------
  for (idx_counter in 1:BENCHMARK_SAMPLES) {
    # Generate data that can be estimated by TMB_GH
    # and is tested on the original parameters from the beginning of this file
    # The goal is to have a dataset that poses no estimation problem,
    # when estimated with guessed initial parameters
    benchmark_model <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_LAMB,
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
    estim_benchmarks_df_lamb <- rbind(estim_benchmarks_df_lamb,
                                      data.frame(time = c(timeDM, timeTMB, timeTMB_G, timeTMB_H, timeTMB_GH),
                                                 m = rep(m, length(PROCEDURES)),
                                                 procedure = PROCEDURES,
                                                 dataset_number = rep(idx_counter, length(PROCEDURES))))
    
    # Benchmark mllk times --------------------------------------------------
    tmb_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
                                  parameters = parvect_benchmark,
                                  MakeADFun_obj = obj_benchmark)
    tmb_g_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
                                    parameters = parvect_benchmark,
                                    MakeADFun_obj = obj_benchmark,
                                    gradient = TRUE)
    tmb_h_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
                                    parameters = parvect_benchmark,
                                    MakeADFun_obj = obj_benchmark,
                                    hessian = TRUE)
    tmb_gh_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
                                     parameters = parvect_benchmark,
                                     MakeADFun_obj = obj_benchmark,
                                     gradient = TRUE,
                                     hessian = TRUE)
    
    param_tmb <- pois.HMM.pn2pw(m = m, lambda = tmb_benchmark$lambda, gamma = tmb_benchmark$gamma)
    param_tmb_g <- pois.HMM.pn2pw(m = m, lambda = tmb_g_benchmark$lambda, gamma = tmb_g_benchmark$gamma)
    param_tmb_h <- pois.HMM.pn2pw(m = m, lambda = tmb_h_benchmark$lambda, gamma = tmb_h_benchmark$gamma)
    param_tmb_gh <- pois.HMM.pn2pw(m = m, lambda = tmb_gh_benchmark$lambda, gamma = tmb_gh_benchmark$gamma)
    model1 <- MakeADFun(TMB_data, param_tmb, DLL = "poi_hmm", silent = TRUE)
    model2 <- MakeADFun(TMB_data, param_tmb_g, DLL = "poi_hmm", silent = TRUE)
    model3 <- MakeADFun(TMB_data, param_tmb_h, DLL = "poi_hmm", silent = TRUE)
    model4 <- MakeADFun(TMB_data, param_tmb_gh, DLL = "poi_hmm", silent = TRUE)
    parvect_mllk <- pois.HMM.pn2pw(m, lambda = tmb_gh_benchmark$lambda, gamma = tmb_gh_benchmark$gamma)
    
    temp <- microbenchmark("DM" = pois.HMM.mllk(parvect_mllk, lamb_data, m),
                           "TMB" = model1$fn(model1$par),
                           "TMB_G" = model2$fn(model2$par),
                           "TMB_H" = model3$fn(model3$par),
                           "TMB_GH" = model4$fn(model4$par),
                           times = 1)
    
    times <- temp$time/10^9
    timeDM <- times[temp$expr == "DM"]
    timeTMB <- times[temp$expr == "TMB"]
    timeTMB_G <- times[temp$expr == "TMB_G"]
    timeTMB_H <- times[temp$expr == "TMB_H"]
    timeTMB_GH <- times[temp$expr == "TMB_GH"]
    mllk_times_df_lamb <- rbind(mllk_times_df_lamb,
                                data.frame(time = c(timeDM, timeTMB, timeTMB_G, timeTMB_H, timeTMB_GH),
                                           m = rep(m, length(PROCEDURES)),
                                           procedure = PROCEDURES,
                                           dataset_number = rep(idx_counter, length(PROCEDURES))))
    
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
    method_comparison_df_lamb <- rbind(method_comparison_df_lamb,
                                       data.frame(time = c(timeBFGS, timeL_BFGS_B, timenlm, timenlminb, timehjn, timemarqLevAlg),
                                                  m = rep(m, length(PROCEDURES_METHOD)),
                                                  procedure = PROCEDURES_METHOD,
                                                  dataset_number = idx_counter))
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
    if (anyNA(ci)) {
      working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
      working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
    } else {
      working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- min(ci)
      working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- max(ci)
    }
  }
  
  # Transform the working parameters into natural ones
  # Lambda (m values)
  conf_int_lamb$Profile.L[which(conf_int_lamb$m == m)][1:m] <- exp(working_conf_int$lower[1:m])
  # Gamma (m^2-m working parameters, m^2 natural ones)
  if (!anyNA(working_conf_int$lower[(m + 1):(m ^ 2)])) {
    natural_gamma <- as.numeric(Gamma_w2n(m, working_conf_int$lower[(m + 1):(m ^ 2)]))
    conf_int_lamb$Profile.L[which(conf_int_lamb$m == m)][(m + 1):(m ^ 2 + m)] <- natural_gamma
  }
  # Lambda (m values)
  conf_int_lamb$Profile.U[which(conf_int_lamb$m == m)][1:m] <- exp(working_conf_int$upper[1:m])
  # Gamma (m^2-m working parameters, m^2 natural ones)
  if (!anyNA(working_conf_int$upper[(m + 1):(m ^ 2)])) {
    natural_gamma <- as.numeric(Gamma_w2n(m, working_conf_int$upper[(m + 1):(m ^ 2)]))
    conf_int_lamb$Profile.U[which(conf_int_lamb$m == m)][(m + 1):(m ^ 2 + m)] <- natural_gamma
  }
  
  # Bootstrap ---------------------------
  bootstrap_lamb <- data.frame()
  lambda <- tmb_CI$lambda
  gamma <- tmb_CI$gamma
  delta <- tmb_CI$delta
  for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
    temp <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_LAMB,
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
    bootstrap_lamb[idx_sample, 1:len_par] <- natural_parameters
  }
  colnames(bootstrap_lamb) <- params_names_latex
  q <- apply(bootstrap_lamb, 2, quantile.colwise)
  conf_int_lamb$Bootstrap.L[which(conf_int_lamb$m == m)] <- q[1, ]
  conf_int_lamb$Bootstrap.U[which(conf_int_lamb$m == m)] <- q[2, ]
  
  # TMB confidence intervals --------------
  conf_int_lamb$TMB.L[which(conf_int_lamb$m == m)] <- c(tmb_gh$lambda - q95_norm * tmb_gh$lambda_std_error,
                                                        tmb_gh$gamma - q95_norm * tmb_gh$gamma_std_error,
                                                        tmb_gh$delta - q95_norm * tmb_gh$delta_std_error)
  conf_int_lamb$TMB.U[which(conf_int_lamb$m == m)] <- c(tmb_gh$lambda + q95_norm * tmb_gh$lambda_std_error,
                                                        tmb_gh$gamma + q95_norm * tmb_gh$gamma_std_error,
                                                        tmb_gh$delta + q95_norm * tmb_gh$delta_std_error)
  # Coverage probabilities of the 3 CI methods -----------------
  parameter_names <- paste0(rep("lambda", m), 1:m)
  for (gamma_idx in 1:m ^ 2) {
    # Get row and column indexes for gamma instead of the default
    # columnwise index: the default indexes are 1:m for the 1st column,
    # then (m + 1):(2 * m) for the 2nd, etc...
    row_col_idx <- matrix.col.idx.to.rowcol(gamma_idx, m)
    parameter_names <- c(parameter_names,
                         paste0("gamma", toString(row_col_idx)))
  }
  parameter_names <- c(parameter_names, paste0(rep("delta", m), 1:m))
  coverage_count_profile <- coverage_count_bootstrap <- coverage_count_tmb <- data.frame(parameter = parameter_names,
                                                                                         count = 0,
                                                                                         ratio = 0)
  # Sometimes, there are no profile CIs for gamma because of the need to transform the CIs to their natural form
  # In these cases, we reduce the total
  coverage_adjusted_sample_size_gamma <- COVERAGE_SAMPLES
  
  for (idx_coverage in 1:COVERAGE_SAMPLES) {
    # Generate a data sample where nlminb converges
    # Loop as long as there is an issue with nlminb
    # Estimate a model
    coverage_model <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_LAMB,
                                                         mod = list(m = m,
                                                                    lambda = tmb_CI$lambda,
                                                                    gamma = tmb_CI$gamma,
                                                                    delta = tmb_CI$delta),
                                                         testing_params = list(m = m,
                                                                               lambda = lambda_init,
                                                                               gamma = gamma_init,
                                                                               delta = delta_init),
                                                         std_error = TRUE)
    
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
      if (anyNA(ci)) {
        working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
        working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
      } else {
        working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- min(ci)
        working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- max(ci)
      }
    }
    
    lambda_profile_lower <- exp(working_conf_int$lower[1:m])
    lambda_profile_upper <- exp(working_conf_int$upper[1:m])
    tmp <- as.numeric(Gamma_w2n(m, working_conf_int$lower[(m + 1):(m ^ 2)]))
    tmp2 <- as.numeric(Gamma_w2n(m, working_conf_int$upper[(m + 1):(m ^ 2)]))
    gamma_profile_lower <- pmin(tmp, tmp2)
    gamma_profile_upper <- pmax(tmp, tmp2)
    
    # If the "true" value for lambda is in the CI, then increase the count
    indexes <- which(as.vector(tmb_CI$lambda) >= lambda_profile_lower & as.vector(tmb_CI$lambda) <= lambda_profile_upper)
    coverage_count_profile[indexes, "count"] <- coverage_count_profile[indexes, "count"] + 1
    # Either or both bounds for gamma can be missing due to missing tgamma values
    # If there are missing tgamma values, the confidence interval has 0 or 1 bound
    # and is thus useless for coverage probabilities
    if (anyNA(gamma_profile_lower, recursive = TRUE) || anyNA(gamma_profile_upper, recursive = TRUE)) {
      coverage_adjusted_sample_size_gamma <- coverage_adjusted_sample_size_gamma - 1
    } else {
      # Same as above for gamma
      indexes <- which(as.vector(tmb_CI$gamma) >= gamma_profile_lower & as.vector(tmb_CI$gamma) <= gamma_profile_upper)
      indexes <- indexes + m
      coverage_count_profile[indexes, "count"] <- coverage_count_profile[indexes, "count"] + 1
    }
    
    # Confidence interval bootstrap -----------------------------------
    bootstrap_lamb <- data.frame()
    lambda_coverage <- coverage_model$mod$lambda
    gamma_coverage <- coverage_model$mod$gamma
    delta_coverage <- coverage_model$mod$delta
    for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
      temp <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_HOSP,
                                                 mod = list(m = m,
                                                            lambda = lambda_coverage,
                                                            gamma = gamma_coverage,
                                                            delta = delta_coverage),
                                                 testing_params = list(m = m,
                                                                       lambda = tmb_CI$lambda,
                                                                       gamma = tmb_CI$gamma,
                                                                       delta = tmb_CI$delta))$natural_parameters
      # The values from gamma are taken columnwise
      natural_parameters <- temp
      natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
      bootstrap_lamb[idx_sample, 1:len_par] <- natural_parameters
    }
    # colnames(bootstrap_lamb) <- params_names_latex
    q <- apply(bootstrap_lamb, 2, quantile.colwise)
    # lambda is at positions 1:m
    # gamma is at positions (m + 1):(m + m ^ 2)
    # delta is at positions (m + m ^ 2 + 1):(m ^ 2 + 2 * m)
    indexes <- which(as.vector(tmb_CI$lambda) >= q[1, 1:m] & as.vector(tmb_CI$lambda) <= q[2, 1:m])
    coverage_count_bootstrap[indexes, "count"] <- coverage_count_bootstrap[indexes, "count"] + 1
    indexes <- which(as.vector(tmb_CI$gamma) >= q[1, (m + 1):(m + m ^ 2)] & as.vector(tmb_CI$gamma) <= q[2, (m + 1):(m + m ^ 2)]) + m
    coverage_count_bootstrap[indexes, "count"] <- coverage_count_bootstrap[indexes, "count"] + 1
    indexes <- which(as.vector(tmb_CI$delta) >= q[1, (m + m ^ 2 + 1):(m ^ 2 + 2 * m)] & as.vector(tmb_CI$delta) <= q[2, (m + m ^ 2 + 1):(m ^ 2 + 2 * m)]) + m + m ^ 2
    coverage_count_bootstrap[indexes, "count"] <- coverage_count_bootstrap[indexes, "count"] + 1
    
    # Confidence interval TMB -----------------------------------------
    lambda_tmb_lower <- coverage_model$natural_parameters$lambda - q95_norm * coverage_model$natural_parameters$lambda_std_error
    lambda_tmb_upper <- coverage_model$natural_parameters$lambda + q95_norm * coverage_model$natural_parameters$lambda_std_error
    gamma_tmb_lower <- coverage_model$natural_parameters$gamma - q95_norm * coverage_model$natural_parameters$gamma_std_error
    gamma_tmb_upper <- coverage_model$natural_parameters$gamma + q95_norm * coverage_model$natural_parameters$gamma_std_error
    delta_tmb_lower <- coverage_model$natural_parameters$delta - q95_norm * coverage_model$natural_parameters$delta_std_error
    delta_tmb_upper <- coverage_model$natural_parameters$delta + q95_norm * coverage_model$natural_parameters$delta_std_error
    
    indexes <- which(as.vector(tmb_CI$lambda) >= lambda_tmb_lower & as.vector(tmb_CI$lambda) <= lambda_tmb_upper)
    coverage_count_tmb[indexes, "count"] <- coverage_count_tmb[indexes, "count"] + 1
    indexes <- which(as.vector(tmb_CI$gamma) >= gamma_tmb_lower & as.vector(tmb_CI$gamma) <= gamma_tmb_upper)
    indexes <- indexes + m
    coverage_count_tmb[indexes, "count"] <- coverage_count_tmb[indexes, "count"] + 1
    indexes <- which(as.vector(tmb_CI$delta) >= delta_tmb_lower & as.vector(tmb_CI$delta) <= delta_tmb_upper)
    indexes <- indexes + m + m ^ 2
    coverage_count_tmb[indexes, "count"] <- coverage_count_tmb[indexes, "count"] + 1
    
  }
  coverage_count_profile[1:m, "ratio"] <- coverage_count_profile[1:m, "count"] / COVERAGE_SAMPLES
  coverage_count_profile[(m + 1):(m + m ^ 2), "ratio"] <- coverage_count_profile[(m + 1):(m + m ^ 2), "count"] / coverage_adjusted_sample_size_gamma
  coverage_count_profile[(m + m ^ 2 + 1):(m ^ 2 + 2 * m), "ratio"] <- NA
  
  coverage_count_tmb$ratio <- coverage_count_tmb$count / COVERAGE_SAMPLES
  
  coverage_count_bootstrap$ratio <- coverage_count_bootstrap$count / COVERAGE_SAMPLES
  
  conf_int_lamb[conf_int_lamb$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Profile"] <- coverage_count_profile$ratio
  conf_int_lamb[conf_int_lamb$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Bootstrap"] <- coverage_count_bootstrap$ratio
  conf_int_lamb[conf_int_lamb$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.TMB"] <- coverage_count_tmb$ratio
}

# The profile standard errors for gamma are transformed back to their natural form
# and may not be sorted, so we sort them manually again
for (i in 1:length(conf_int_lamb[, 1])) {
  row <- conf_int_lamb[i, c("Profile.L", "Profile.U")]
  conf_int_lamb[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
}
mllk_values_lamb$m <- as.factor(as.numeric(mllk_values_lamb$m))
mllk_values_lamb$mllk <- as.numeric(mllk_values_lamb$mllk)
mllk_values_lamb$AIC <- as.numeric(mllk_values_lamb$AIC)
mllk_values_lamb$BIC <- as.numeric(mllk_values_lamb$BIC)
conf_int_lamb$m <- as.integer(conf_int_lamb$m)

# Remove gamma and delta when m = 1 because they're pointless
if (1 %in% M_LIST_LAMB) {
  conf_int_lamb <- conf_int_lamb[- c(2,3), ]
}

estim_benchmarks_df_lamb$m <- factor(estim_benchmarks_df_lamb$m, levels = M_LIST_LAMB)