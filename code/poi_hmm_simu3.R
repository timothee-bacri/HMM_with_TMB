# Prepare the data & parameters, then estimate for different numbers of hidden states
for (idx in 1:length(M_LIST_SIMU3)) {
  set.seed(18)
  # Parameters and covariates --------------------------
  m <- M_LIST_SIMU3[idx]
  if (m == 2) {
    true_gamma <- matrix(c(0.95, 0.05,
                           0.15, 0.85),
                         byrow = TRUE,
                         nrow = m,
                         ncol = m)
  } else if (m == 3) {
    true_gamma <- matrix(c(0.950, 0.025, 0.025,
                           0.050, 0.900, 0.050,
                           0.075, 0.075, 0.850),
                         byrow = TRUE,
                         nrow = m,
                         ncol = m)
  } else if (m == 4) {
    true_gamma <- matrix(c(0.85, 0.05, 0.05, 0.05,
                           0.05, 0.85, 0.05, 0.05,
                           0.05, 0.10, 0.80, 0.05,
                           0.25, 0.25, 0.25, 0.25),
                           # 0.05, 0.05, 0.05, 0.85),
                         byrow = TRUE,
                         nrow = m,
                         ncol = m)
  }
  true_lambda <- seq(1,
                     16,
                     length.out = m)
  true_delta <- stat.dist(true_gamma)
  
  simu3_data <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                   mod = list(m = m,
                                                              lambda = true_lambda,
                                                              gamma = true_gamma,
                                                              delta = true_delta),
                                                   testing_params = list(m = m,
                                                                         lambda = true_lambda,
                                                                         gamma = true_gamma,
                                                                         delta = true_delta))$data
  
  
  # Parameters & covariates for TMB ------------------
  working_true_params <- pois.HMM.pn2pw(m = m,
                                        lambda = true_lambda,
                                        gamma = true_gamma,
                                        delta = true_delta)
  TMB_data <- list(x = simu3_data,
                   m = m)
  
  # Estimation ------------------------------------
  dm <- DM.estimate(x = simu3_data,
                    m = m,
                    lambda0 = true_lambda,
                    gamma0 = true_gamma,
                    delta0 = true_delta)
  tmb <- TMB.estimate(TMB_data = TMB_data,
                      parameters = working_true_params)
  tmb_g <- TMB.estimate(TMB_data = TMB_data,
                        parameters = working_true_params,
                        gradient = TRUE)
  tmb_h <- TMB.estimate(TMB_data = TMB_data,
                        parameters = working_true_params,
                        hessian = TRUE)
  tmb_gh <- TMB.estimate(TMB_data = TMB_data,
                         parameters = working_true_params,
                         gradient = TRUE,
                         hessian = TRUE,
                         std_error = TRUE)
  
  # If one doesn't converge successfully, stop
  # if (dm$convergence != 0) {
  #   stop(paste("dm didn't converge properly, simu3 dataset, m =",
  #              m))
  # }
  if (tmb$convergence != 0) {
    stop(paste("tmb didn't converge properly, simu3 dataset, m =",
               m))
  }
  if (tmb_g$convergence != 0) {
    stop(paste("tmb_g didn't converge properly, simu3 dataset, m =",
               m))
  }
  if (tmb_h$convergence != 0) {
    stop(paste("tmb_h didn't converge properly, simu3 dataset, m =",
               m))
  }
  if (tmb_gh$convergence != 0) {
    stop(paste("tmb_gh didn't converge properly, simu3 dataset, m =",
               m))
  }
  
  # Creating variables for the CIs -----------------
  tmb_CI <- tmb_gh
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
  indices <- (length(conf_int_simu3$m) + 1):(length(conf_int_simu3$m) + len_par)
  conf_int_simu3[indices, "m"] <- m
  conf_int_simu3[indices, "Parameter"] <- params_names_latex
  # Reminder, PARAMS_NAMES contains c("lambda", "gamma", "delta")
  conf_int_simu3[indices, "Estimate"] <- unlist(tmb_CI[PARAMS_NAMES])
  conf_int_simu3[indices, "True.value"] <- as.numeric(c(true_lambda,
                                                        true_gamma,
                                                        true_delta))
  
  param_tmb_CI <- pois.HMM.pn2pw(m = m,
                                 lambda = tmb_CI$lambda,
                                 gamma = tmb_CI$gamma)
  
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
  
  lambda_indices <- 1:m
  gamma_indices <- max(lambda_indices) + 1:(m ^ 2)
  delta_indices <- max(gamma_indices) + 1:m
  tgamma_indices <- (m + 1):(m ^ 2)
  
  # Benchmarks ------------
  set.seed(19)
  if (BENCHMARK_SAMPLES != 0) {
    sms("Simu3 Benchmark start")
    for (idx_counter in 1:BENCHMARK_SAMPLES) {
      # Generate data that can be estimated by TMB_GH
      # and is tested on the slightly off parameters from the beginning of this file
      # The goal is to have a dataset that poses no estimation problem,
      # when estimated with guessed initial parameters
      benchmark_model <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                            mod = list(m = m,
                                                                       lambda = tmb_CI$lambda,
                                                                       gamma = tmb_CI$gamma,
                                                                       delta = tmb_CI$delta),
                                                            testing_params = list(m = m,
                                                                                  lambda = true_lambda,
                                                                                  gamma = true_gamma,
                                                                                  delta = true_delta),
                                                            test_marqLevAlg = TRUE)
      benchmark_data <- benchmark_model$data
      # Benchmark all different combinations of gradient and hessians with DM ----------------
      # Parameters & covariates for DM and TMB
      TMB_benchmark_data <- list(x = benchmark_data,
                                 m = m)
      obj_benchmark <- MakeADFun(TMB_benchmark_data,
                                 working_true_params,
                                 DLL = "poi_hmm",
                                 silent = TRUE)
      parvect_benchmark_TMB <- pois.HMM.pn2pw(m = m,
                                          lambda = true_lambda,
                                          gamma = true_gamma,
                                          delta = true_delta)
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
                                                                 working_true_params,
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
      
      estim_benchmarks_df_simu3 <- rbind(estim_benchmarks_df_simu3,
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
      
      # Benchmark mllk times --------------------------------------------------
      tmb_gh_benchmark <- TMB.estimate(TMB_data = TMB_benchmark_data,
                                       parameters = parvect_benchmark_TMB,
                                       gradient = TRUE,
                                       hessian = TRUE)
      
      param_tmb_gh <- pois.HMM.pn2pw(m = m,
                                     lambda = tmb_gh_benchmark$lambda,
                                     gamma = tmb_gh_benchmark$gamma)
      model4 <- MakeADFun(TMB_data,
                          param_tmb_gh,
                          DLL = "poi_hmm",
                          silent = TRUE)
      parvect_mllk <- pois.HMM.pn2pw(m,
                                     lambda = tmb_gh_benchmark$lambda,
                                     gamma = tmb_gh_benchmark$gamma)
      
      temp <- microbenchmark("DM" = pois.HMM.mllk(parvect_mllk,
                                                  simu3_data,
                                                  m),
                             "TMB_GH" = model4$fn(model4$par),
                             times = 1)
      
      times <- temp$time / 10^9
      timeDM <- times[temp$expr == "DM"]
      timeTMB_GH <- times[temp$expr == "TMB_GH"]
      mllk_times_df_simu3 <- rbind(mllk_times_df_simu3,
                                   data.frame(time = c(timeDM,
                                                       timeTMB_GH),
                                              m = rep(m,
                                                      2),
                                              procedure = PROCEDURES[c(1, 5)],
                                              dataset_number = rep(idx_counter,
                                                                   2)))
      
      # Benchmark different optimization methods ----------------------------------------------
      temp <- microbenchmark("BFGS" = optim(par = obj_benchmark$par,
                                            fn = obj_benchmark$fn,
                                            gr = obj_benchmark$gr,
                                            method = "BFGS",
                                            control = ctrl)$convergence==0,
                             "L-BFGS-B" = optim(par = obj_benchmark$par,
                                                fn = obj_benchmark$fn,
                                                gr = obj_benchmark$gr,
                                                method = "L-BFGS-B",
                                                control = ctrl)$convergence==0,
                             "nlm" = nlm(f = nlmfn,
                                         p = obj_benchmark$par,
                                         obj_benchmark,
                                         iterlim = 10000)$code==1,
                             "nlminb" = nlminb(start = obj_benchmark$par,
                                               objective = obj_benchmark$fn,
                                               gradient = obj_benchmark$gr,
                                               hessian = obj_benchmark$he)$convergence==0,
                             "hjn" = hjn(par = obj_benchmark$par,
                                         fn = obj_benchmark$fn,
                                         lower = HJN_LOWER_LIMIT,
                                         upper = HJN_UPPER_LIMIT)$convergence==0,
                             "marqLevAlg" = marqLevAlg(b = obj_benchmark$par,
                                                       fn = obj_benchmark$fn,
                                                       gr = obj_benchmark$gr,
                                                       hess = obj_benchmark$he,
                                                       maxiter = 10000)$istop==1,
                             times = 1,
                             check = "equal",
                             setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data,
                                                                 working_true_params,
                                                                 DLL = "poi_hmm",
                                                                 silent = TRUE))
      
      times <- temp$time / 10^9
      timeBFGS <- times[temp$expr == "BFGS"]
      timeL_BFGS_B <- times[temp$expr == "L-BFGS-B"]
      timenlm <- times[temp$expr == "nlm"]
      timenlminb <- times[temp$expr == "nlminb"]
      timehjn <- times[temp$expr == "hjn"]
      timemarqLevAlg <- times[temp$expr == "marqLevAlg"]
      method_comparison_df_simu3 <- rbind(method_comparison_df_simu3,
                                          data.frame(time = c(timeBFGS,
                                                              timeL_BFGS_B,
                                                              timenlm,
                                                              timenlminb,
                                                              timehjn,
                                                              timemarqLevAlg),
                                                     m = rep(m,
                                                             length(PROCEDURES_METHOD)),
                                                     procedure = PROCEDURES_METHOD,
                                                     dataset_number = idx_counter))
    }
    sms("Simu3 Benchmark end")
  }
  
  # Profiling the likelihood --------------------------
  
  # SEQUENTIAL/NORMAL
  # working_conf_int <- data.frame(w_parameter = w_params_names,
  #                                lower = rep(NA,
  #                                            len_w_par),
  #                                upper = rep(NA,
  #                                            len_w_par),
  #                                stringsAsFactors = FALSE)
  # for (idx_param in 1:len_w_par) {
  #   profile <- tmbprofile(obj = tmb_CI$obj,
  #                         name = idx_param,
  #                         trace = FALSE)
  #   
  #   ci <- tryCatch({
  #     confint(profile)
  #   },
  #   error = function(e){
  #     return(rep(NA,
  #                2))
  #   })
  #   
  #   w_param_name <- w_params_names[idx_param]
  #   working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
  #   working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
  # }
  # 
  
  # PARALLEL
  registerDoParallel(cores = CORES)
  working_conf_int <- foreach (idx_param = 1:len_w_par,
                               .packages = "TMB",
                               .inorder = TRUE,
                               .combine = rbind) %dopar% {
                                 # TMB::compile("code/poi_hmm.cpp")
                                 dyn.load(dynlib("code/poi_hmm"))
                                 profile <- tmbprofile(obj = tmb_CI$obj,
                                                       name = idx_param,
                                                       trace = FALSE)
                                 
                                 ci <- tryCatch({
                                   confint(profile)
                                 },
                                 error = function(e){
                                   return(rep(NA,
                                              2))
                                 })
                                 
                                 # w_param_name <- w_params_names[idx_param]
                                 # working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
                                 # working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
                                 dyn.unload(dynlib("code/poi_hmm"))
                                 return(ci)
                               }
  stopImplicitCluster()
  working_conf_int <- as.data.frame(working_conf_int)
  rownames(working_conf_int) <- w_params_names
  
  # Transform the working parameters into natural ones
  # Lambda (m values)
  conf_int_simu3$Profile.L[which(conf_int_simu3$m == m)][lambda_indices] <- exp(working_conf_int$lower[lambda_indices])
  # Gamma (m^2-m working parameters, m^2 natural ones)
  if (!anyNA(working_conf_int$lower[tgamma_indices])) {
    natural_gamma <- as.numeric(gamma.w2n(m,
                                          working_conf_int$lower[tgamma_indices]))
    conf_int_simu3$Profile.L[which(conf_int_simu3$m == m)][gamma_indices] <- natural_gamma
  }
  # Lambda (m values)
  conf_int_simu3$Profile.U[which(conf_int_simu3$m == m)][lambda_indices] <- exp(working_conf_int$upper[lambda_indices])
  # Gamma (m^2-m working parameters, m^2 natural ones)
  if (!anyNA(working_conf_int$upper[tgamma_indices])) {
    natural_gamma <- as.numeric(gamma.w2n(m,
                                          working_conf_int$upper[tgamma_indices]))
    conf_int_simu3$Profile.U[which(conf_int_simu3$m == m)][gamma_indices] <- natural_gamma
  }
  
  # Bootstrap ---------------------------
  set.seed(20)
  if (BOOTSTRAP_SAMPLES != 0) {
    sms("Simu3 Bootstrap start")
    registerDoParallel(cores = CORES)
    bootstrap_simu3 <- foreach (idx_sample = 1:BOOTSTRAP_SAMPLES,
                                .packages = "TMB",
                                .combine = rbind) %dopar% {
                                  # TMB::compile("code/poi_hmm.cpp")
                                  dyn.load(dynlib("code/poi_hmm"))
                                  temp <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                                             mod = list(m = m,
                                                                                        lambda = tmb_CI$lambda,
                                                                                        gamma = tmb_CI$gamma,
                                                                                        delta = tmb_CI$delta),
                                                                             testing_params = list(m = m,
                                                                                                   lambda = true_lambda,
                                                                                                   gamma = true_gamma,
                                                                                                   delta = true_delta))$natural_parameters
                                  # The values from gamma are taken columnwise
                                  natural_parameters <- temp
                                  natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
                                  dyn.unload(dynlib("code/poi_hmm"))
                                  return(natural_parameters)
                                }
    stopImplicitCluster()
    colnames(bootstrap_simu3) <- params_names_latex
    q <- apply(bootstrap_simu3,
               2,
               quantile.colwise)
    conf_int_simu3$Bootstrap.L[which(conf_int_simu3$m == m)] <- q[1, ]
    conf_int_simu3$Bootstrap.U[which(conf_int_simu3$m == m)] <- q[2, ]
    sms("Simu3 Bootstrap end")
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
  conf_int_simu3$TMB.L[which(conf_int_simu3$m == m)] <- c(lambda_L,
                                                          gamma_L,
                                                          delta_L)
  # no upper bound on lambda
  # gamma must be 1 or less
  gamma_U <- pmin(1,
                  tmb_gh$gamma + q95_norm * tmb_gh$gamma_std_error)
  # delta must be 1 or less
  delta_U <- pmin(1,
                  tmb_gh$delta + q95_norm * tmb_gh$delta_std_error)
  conf_int_simu3$TMB.U[which(conf_int_simu3$m == m)] <- c(tmb_gh$lambda + q95_norm * tmb_gh$lambda_std_error,
                                                          gamma_U,
                                                          delta_U)
  # Coverage probabilities of the 3 CI methods -----------------
  set.seed(21)
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
  coverage_count_profile <- coverage_count_bootstrap <- coverage_count_tmb <- data.frame(parameter = parameter_names,
                                                                                         count = 0,
                                                                                         ratio = 0)
  idx_coverage <- 0
  sms("Simu3 Coverage start")
  begin_coverage <- Sys.time()
  while (idx_coverage < COVERAGE_SAMPLES) {
    idx_coverage <- idx_coverage + 1
    # Generate a data sample where nlminb converges
    # Loop as long as there is an issue with nlminb
    # Estimate a model
    # Unlike with the other datasets, we know the true parameters of this one
    coverage_model <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                         mod = list(m = m,
                                                                    lambda = true_lambda,
                                                                    gamma = true_gamma,
                                                                    delta = true_delta),
                                                         testing_params = list(m = m,
                                                                               lambda = true_lambda,
                                                                               gamma = true_gamma,
                                                                               delta = true_delta),
                                                         std_error = TRUE)
    
    # Save the occurrences of failures to generate a sample for which parameters can be estimated
    for (reason in c("state_number", "TMB_null", "TMB_converge", "TMB_G_null",
                     "TMB_G_converge", "TMB_H_null", "TMB_H_converge", "TMG_GH_null",
                     "TMG_GH_converge", "marqLevAlg_converge", "NA_value")) {
      coverage_skips_simu3[coverage_skips_simu3$m == m, reason] <- coverage_skips_simu3[coverage_skips_simu3$m == m, reason] + coverage_model$failure[reason]
    }
    
    # Confidence interval profiling -------------------------------------
    registerDoParallel(cores = CORES)
    working_conf_int <- foreach (idx_param = 1:len_w_par,
                                 .packages = "TMB",
                                 .inorder = TRUE,
                                 .combine = rbind) %dopar% {
                                   # TMB::compile("code/poi_hmm.cpp")
                                   dyn.load(dynlib("code/poi_hmm"))
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
                                   dyn.unload(dynlib("code/poi_hmm"))
                                   return(ci)
                                 }
    stopImplicitCluster()
    working_conf_int <- as.data.frame(working_conf_int)
    rownames(working_conf_int) <- w_params_names
    
    lambda_profile_lower <- exp(working_conf_int$lower[lambda_indices])
    lambda_profile_upper <- exp(working_conf_int$upper[lambda_indices])
    gamma_profile_lower <- as.numeric(gamma.w2n(m,
                                                working_conf_int$lower[tgamma_indices]))
    gamma_profile_upper <- as.numeric(gamma.w2n(m,
                                                working_conf_int$upper[tgamma_indices]))
    
    # If profiling doesn't yield results for all parameters, try a new coverage sample
    if (anyNA(c(lambda_profile_lower,
                lambda_profile_upper,
                gamma_profile_lower,
                gamma_profile_upper),
              recursive = TRUE)) {
      idx_coverage <- idx_coverage - 1
      coverage_skips_simu3[coverage_skips_simu3$m == m, "profile"] <- coverage_skips_simu3[coverage_skips_simu3$m == m, "profile"] + 1
      next
    }
    # If the true value for lambda is in the CI, then increase the count
    real_lambda_profile_lower <- pmin(lambda_profile_lower,
                                      lambda_profile_upper)
    real_lambda_profile_upper <- pmax(lambda_profile_lower,
                                      lambda_profile_upper)
    indices <- which(true_lambda >= real_lambda_profile_lower & true_lambda <= real_lambda_profile_upper)
    coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
    
    # Same for gamma
    real_gamma_profile_lower <- pmin(gamma_profile_lower,
                                     gamma_profile_upper)
    real_gamma_profile_upper <- pmax(gamma_profile_lower,
                                     gamma_profile_upper)
    indices <- which(as.vector(true_gamma) >= real_gamma_profile_lower & as.vector(true_gamma) <= real_gamma_profile_upper)
    indices <- indices + m
    coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
    
    # Confidence interval bootstrap -----------------------------------
    if (BOOTSTRAP_SAMPLES != 0) {
      registerDoParallel(cores = CORES)
      bootstrap_simu3 <- foreach (idx_sample = 1:BOOTSTRAP_SAMPLES,
                                  .packages = "TMB",
                                  .combine = rbind) %dopar% {
                                    # TMB::compile("code/poi_hmm.cpp")
                                    dyn.load(dynlib("code/poi_hmm"))
                                    temp <- pois.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                                               mod = list(m = m,
                                                                                          lambda = coverage_model$natural_parameters$lambda,
                                                                                          gamma = coverage_model$natural_parameters$gamma,
                                                                                          delta = coverage_model$natural_parameters$delta),
                                                                               testing_params = list(m = m,
                                                                                                     lambda = true_lambda,
                                                                                                     gamma = true_gamma,
                                                                                                     delta = true_delta))$natural_parameters
                                    # The values from gamma are taken columnwise
                                    natural_parameters <- temp
                                    natural_parameters <- unlist(natural_parameters[PARAMS_NAMES])
                                    dyn.unload(dynlib("code/poi_hmm"))
                                    return(natural_parameters)
                                  }
      stopImplicitCluster()
      q <- apply(bootstrap_simu3,
                 2,
                 quantile.colwise)
      indices <- which(as.vector(true_lambda) >= q[1, lambda_indices] & as.vector(true_lambda) <= q[2, lambda_indices])
      coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
      indices <- which(as.vector(true_gamma) >= q[1, gamma_indices] & as.vector(true_gamma) <= q[2, gamma_indices]) + m
      coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
      indices <- which(as.vector(true_delta) >= q[1, delta_indices] & as.vector(true_delta) <= q[2, delta_indices]) + m + m ^ 2
      coverage_count_bootstrap[indices, "count"] <- coverage_count_bootstrap[indices, "count"] + 1
    }
    
    # Confidence interval TMB -----------------------------------------
    lambda_tmb_lower <- coverage_model$natural_parameters$lambda - q95_norm * coverage_model$natural_parameters$lambda_std_error
    lambda_tmb_upper <- coverage_model$natural_parameters$lambda + q95_norm * coverage_model$natural_parameters$lambda_std_error
    gamma_tmb_lower <- coverage_model$natural_parameters$gamma - q95_norm * coverage_model$natural_parameters$gamma_std_error
    gamma_tmb_upper <- coverage_model$natural_parameters$gamma + q95_norm * coverage_model$natural_parameters$gamma_std_error
    delta_tmb_lower <- coverage_model$natural_parameters$delta - q95_norm * coverage_model$natural_parameters$delta_std_error
    delta_tmb_upper <- coverage_model$natural_parameters$delta + q95_norm * coverage_model$natural_parameters$delta_std_error
    
    indices <- which(as.vector(true_lambda) >= lambda_tmb_lower & as.vector(true_lambda) <= lambda_tmb_upper)
    coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
    indices <- which(as.vector(true_gamma) >= gamma_tmb_lower & as.vector(true_gamma) <= gamma_tmb_upper)
    indices <- indices + m
    coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
    indices <- which(as.vector(true_delta) >= delta_tmb_lower & as.vector(true_delta) <= delta_tmb_upper)
    indices <- indices + m + m ^ 2
    coverage_count_tmb[indices, "count"] <- coverage_count_tmb[indices, "count"] + 1
    
    percent <- idx_coverage / COVERAGE_SAMPLES * 100
    if (percent %% 20 == 0 && percent != 100) {
      now <- Sys.time()
      percentage <- idx_coverage / COVERAGE_SAMPLES
      duration <- now - begin_coverage
      time_to_completion <- duration / percentage
      
      expected_end_time <- begin_coverage + time_to_completion
      
      sms(paste0("Simu3 ", percent,
                 "%\neta: ", round(expected_end_time - now, 1), " ", units(expected_end_time - now),
                 "\nend at: ", expected_end_time))
    }
  }
  sms("Simu3 Coverage end")
  coverage_count_profile[lambda_indices, "ratio"] <- coverage_count_profile[lambda_indices, "count"] / COVERAGE_SAMPLES
  coverage_count_profile[gamma_indices, "ratio"] <- coverage_count_profile[gamma_indices, "count"] / COVERAGE_SAMPLES
  coverage_count_profile[delta_indices, "ratio"] <- NA # delta is not a parameter for us, so it has no profile CI
  
  coverage_count_tmb$ratio <- coverage_count_tmb$count / COVERAGE_SAMPLES
  
  coverage_count_bootstrap$ratio <- coverage_count_bootstrap$count / COVERAGE_SAMPLES
  
  conf_int_simu3[conf_int_simu3$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Profile"] <- coverage_count_profile$ratio * 100
  conf_int_simu3[conf_int_simu3$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Bootstrap"] <- coverage_count_bootstrap$ratio * 100
  conf_int_simu3[conf_int_simu3$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.TMB"] <- coverage_count_tmb$ratio * 100
  
  # Fixes -------------------------
  # Fix label switching in conf_int_simu3
  ordered_params <- pois.HMM.label.order(m = m,
                                         lambda = true_lambda,
                                         gamma = true_gamma,
                                         delta = true_delta)
  
  new_lambda_indices <- ordered_params$ordered_lambda_indices
  new_gamma_indices <- ordered_params$ordered_gamma_vector_indices
  new_delta_indices <- ordered_params$ordered_delta_indices
  
  conf_int_simu3[conf_int_simu3$m == m, - 2][lambda_indices, ] <- conf_int_simu3[conf_int_simu3$m == m, - 2][lambda_indices, ][new_lambda_indices, ]
  conf_int_simu3[conf_int_simu3$m == m, - 2][gamma_indices, ] <- conf_int_simu3[conf_int_simu3$m == m, - 2][gamma_indices, ][new_gamma_indices, ]
  conf_int_simu3[conf_int_simu3$m == m, - 2][delta_indices, ] <- conf_int_simu3[conf_int_simu3$m == m, - 2][delta_indices, ][new_delta_indices, ]
  
  # Reorder the TPM row-wise instead of column-wise
  # Lexicographical parameter sort for gamma (sort on the parameter name)
  new_gamma_indices_truncated_table <- order(conf_int_simu3[conf_int_simu3$m == m, ][gamma_indices, "Parameter"])
  # Replace rows by sorted rows
  conf_int_simu3[conf_int_simu3$m == m, ][gamma_indices, ] <- conf_int_simu3[conf_int_simu3$m == m, ][gamma_indices, ][new_gamma_indices_truncated_table, ]
}

# The profile CIs may not be sorted, so we sort them manually
for (i in 1:length(conf_int_simu3[, 1])) {
  row <- conf_int_simu3[i, c("Profile.L", "Profile.U")]
  conf_int_simu3[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
}
conf_int_simu3$m <- as.integer(conf_int_simu3$m)

estim_benchmarks_df_simu3$m <- factor(estim_benchmarks_df_simu3$m,
                                      levels = M_LIST_SIMU3)
