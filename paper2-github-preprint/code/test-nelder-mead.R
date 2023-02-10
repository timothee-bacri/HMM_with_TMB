# https://github.com/mclements/rstpm2/blob/master/src/c_optim.cpp#L37
# int trace = 0
# int maxit = 500
# double abstol = -Inf
# double reltol = sqrt(.Machine$double.eps) = 1.490116e-08
# double alpha = 1
# double beta = 0.5
# double gamma = 2
# double epshess = 
# bool hessianp = 

source("code/setup_parameters.R")
source("code/packages.R")
sp500_scaled_data <- 100 * sp500_data
library(tictoc)
library(progressr)
library(TMB)
library(data.table)
dyn.load(dynlib("code/poi_hmm"))
dyn.load(dynlib("code/norm_hmm"))
it <- 1:400
SEED <- 1
BENCHMARK_NELDER_MEAD <- 100
ITERATIONS_BENCHMARK_HYBRID <- 1:200
SAMPLE_SIZE_INITIAL_CONDITIONS <- 1000
# BENCHMARK_NELDER_MEAD <- 3
# ITERATIONS_BENCHMARK_HYBRID <- 1:3
# SAMPLE_SIZE_INITIAL_CONDITIONS <- 20
handlers(handler_progress(
  format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed. ETA: :eta",
  # width    = 60,
  complete = "+"
))

# Graphical tests iteration number threshold for neldermead alone
if (FALSE) {
  SEED <- 1
  # Prepare the data & parameters, then estimate for different numbers of hidden states
  set.seed(SEED)
  # Parameters and covariates tinn --------------------------
  m <- M_LIST_TINN
  gamma <- matrix(0.4 / (m - 1),
                  nrow = m,
                  ncol = m)
  diag(gamma) <- 0.6
  # lambda <- seq(quantile(tinn_data,
  #                             0.1),
  #                    quantile(tinn_data,
  #                             0.9),
  #                    length.out = m)
  lambda <- c(0.5, 0.6)
  delta <- stat.dist(gamma)
  
  # Parameters & covariates for TMB ------------------
  working_params <- pois.HMM.pn2pw(m,
                                   lambda,
                                   gamma)
  
  TMB_data <- list(x = tinn_data,
                   m = m)
  
  # Iterations -------------------
  tinn_nll <- TMB.estimate(TMB_data = TMB_data,
                           working_parameters = working_params,
                           gradient = TRUE,
                           hessian = TRUE)$nll
  
  obj <- MakeADFun(TMB_data,
                   working_params,
                   DLL = "poi_hmm",
                   silent = TRUE)
  nll <- conv <- c()
  for (i in it) {
    neldermead <- optim(par = obj$par,
                        fn = obj$fn,
                        method = "Nelder-Mead",
                        control = list(maxit = i)) # 500
    # res <- pois.HMM.pw2pn(m, neldermead$par)
    nll[i] <- neldermead$value
    conv[i] <- neldermead$convergence == 0
  }
  results <- data.frame(nll = nll, conv = conv, it = it, true_nll = tinn_nll)
  ggplot(results, aes(x = it, y = nll)) +
    geom_point(aes(fill = conv, color = conv)) +
    geom_line(aes(color = conv)) +
    geom_line(aes(y = true_nll)) +
    ggtitle("TYT")
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Prepare the data & parameters, then estimate for different numbers of hidden states
  set.seed(SEED)
  # Parameters and covariates simu3 --------------------------
  m <- M_LIST_SIMU3
  if (m == 2) {
    simu3_true_gamma <- matrix(c(0.8, 0.2,
                                 0.25, 0.75),
                               byrow = TRUE,
                               nrow = m,
                               ncol = m)
  } else if (m == 3) {
    simu3_true_gamma <- matrix(c(0.950, 0.025, 0.025,
                                 0.050, 0.900, 0.050,
                                 0.075, 0.075, 0.850),
                               byrow = TRUE,
                               nrow = m,
                               ncol = m)
  } else if (m == 4) {
    simu3_true_gamma <- matrix(c(0.850, 0.050, 0.050, 0.05,
                                 0.050, 0.850, 0.050, 0.05,
                                 0.050, 0.100, 0.800, 0.05,
                                 0.034, 0.033, 0.033, 0.90),
                               byrow = TRUE,
                               nrow = m,
                               ncol = m)
  }
  simu3_true_mu <- seq(-5,
                       5,
                       length.out = m)
  simu3_true_sigma <- seq(0.5, 1, length.out = m)
  simu3_true_delta <- stat.dist(simu3_true_gamma)
  
  simu3_data <- norm.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3+200,
                                                   mod = list(m = m,
                                                              mu = simu3_true_mu,
                                                              sigma = simu3_true_sigma,
                                                              gamma = simu3_true_gamma,
                                                              delta = simu3_true_delta),
                                                   debug_message = "simu3_data")
  simu3_states <- simu3_data$states
  simu3_data <- simu3_data$data
  
  simu3_nll <- norm.HMM.mllk(parvect = norm.HMM.pn2pw(m = m,
                                                      mu = simu3_true_mu,
                                                      sigma = simu3_true_sigma,
                                                      gamma = simu3_true_gamma),
                             x_alias = simu3_data,
                             m_alias = m)
  
  # Different parameters from the true ones, to begin estimation
  gamma_init <- matrix(0.7 / (m - 1),
                       nrow = m,
                       ncol = m)
  diag(gamma_init) <- 0.8
  mu_init <- seq(-2,
                 2,
                 length.out = m)
  sigma_init <- seq(0.2, 0.3, length.out = m)
  delta_init <- stat.dist(gamma_init)
  
  # Parameters & covariates for TMB ------------------
  working_params_init <- norm.HMM.pn2pw(m,
                                        mu_init,
                                        sigma_init,
                                        gamma_init)
  # working_params_init$true_lambda <- true_lambda
  # working_params_init$simu3_true_gamma <- simu3_true_gamma
  # working_params_init$simu3_true_delta <- simu3_true_delta
  
  TMB_data <- list(x = simu3_data,
                   m = m)
  
  # Iterations -------------------
  obj <- MakeADFun(TMB_data,
                   working_params_init,
                   DLL = "norm_hmm",
                   silent = TRUE)
  nll <- conv <- c()
  for (i in it) {
    neldermead <- optim(par = obj$par,
                        fn = obj$fn,
                        method = "Nelder-Mead",
                        control = list(maxit = i)) # 500
    # res <- norm.HMM.pw2pn(m, neldermead$par)
    nll[i] <- neldermead$value
    conv[i] <- neldermead$convergence == 0
  }
  results <- data.frame(nll = nll, conv = conv, it = it, true_nll = simu3_nll)
  ggplot(results, aes(x = it, y = nll)) +
    geom_point(aes(fill = conv, color = conv)) +
    geom_line(aes(color = conv)) +
    geom_line(aes(y = simu3_nll)) +
    ggtitle(paste("Simu norm", length(simu3_data), "data"))
  
  
  
  
  
  
  
  
  
  
  
  
  # Prepare the data & parameters, then estimate for different numbers of hidden states
  set.seed(SEED)
  # Parameters and covariates sp500 --------------------------
  m <- 2
  gamma <- matrix(0.4 / (m - 1),
                  nrow = m,
                  ncol = m)
  diag(gamma) <- 0.6
  # lambda <- seq(quantile(tinn_data,
  #                             0.1),
  #                    quantile(tinn_data,
  #                             0.9),
  #                    length.out = m)
  mu <- c(0.5, 1.5) * mean(sp500_scaled_data)
  sigma <- c(0.5, 1.5) * sd(sp500_scaled_data)
  delta <- stat.dist(gamma)
  
  # Parameters & covariates for TMB ------------------
  working_params <- norm.HMM.pn2pw(m,
                                   mu,
                                   sigma,
                                   gamma)
  
  TMB_data <- list(x = sp500_scaled_data,
                   m = m)
  
  # Iterations -------------------
  sp500_nll <- norm.TMB.estimate(TMB_data = TMB_data,
                                 working_parameters = working_params,
                                 gradient = TRUE,
                                 hessian = TRUE)$nll
  
  obj <- MakeADFun(TMB_data,
                   working_params,
                   DLL = "norm_hmm",
                   silent = TRUE)
  nll <- conv <- c()
  for (i in it) {
    neldermead <- optim(par = obj$par,
                        fn = obj$fn,
                        method = "Nelder-Mead",
                        control = list(maxit = i)) # 500
    # res <- pois.HMM.pw2pn(m, neldermead$par)
    nll[i] <- neldermead$value
    conv[i] <- neldermead$convergence == 0
  }
  results <- data.frame(nll = nll, conv = conv, it = it, true_nll = sp500_nll)
  ggplot(results, aes(x = it, y = nll)) +
    geom_point(aes(fill = conv, color = conv)) +
    geom_line(aes(color = conv)) +
    geom_line(aes(y = true_nll)) +
    ggtitle("SP500")
  
  
}

# Time and compare nll of neldermead alone & hybrid on SP500
if (FALSE) {
  with_progress({
    time_iterations_table <- data.frame(iterations = integer(),
                                        data = character(),
                                        method = character(),
                                        time = character())
    p <- progressor(along = 1:(4 * max(ITERATIONS_BENCHMARK_HYBRID)))
    for (hybrid_iterations in ITERATIONS_BENCHMARK_HYBRID) {
      # Prepare the data & parameters, then estimate for different numbers of hidden states
      # tinn --------------------------
      if (TRUE) {
        set.seed(SEED)
        CONTROL_ARGS$`Nelder-Mead`$maxit <- 10000
        # Parameters and covariates tinn --------------------------
        m <- M_LIST_TINN
        DATA <- tinn_data
        lambda <- seq(quantile(DATA, 0.1, na.rm = TRUE),
                      quantile(DATA, 0.9, na.rm = TRUE),
                      length.out = m)
        gamma <- matrix(0.2 / (m - 1),
                        nrow = m,
                        ncol = m)
        diag(gamma) <- 0.8
        delta <- stat.dist(gamma)
        
        # Parameters & covariates for TMB ------------------
        working_params <- pois.HMM.pn2pw(m = m,
                                         lambda = lambda,
                                         gamma = gamma,
                                         delta = delta)
        TMB_data <- list(x = DATA,
                         m = m)
        
        # Benchmark different optimization methods (ensure high enough iterations to get convergence) -----------------------------
        temp_microbenchmark <- microbenchmark("Nelder-Mead" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "poi_hmm",
                           silent = TRUE)
          optim(par = obj$par,
                fn = obj$fn,
                method = "Nelder-Mead",
                control = CONTROL_ARGS$`Nelder-Mead`)
        },
        
        "nlminb" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "poi_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        
        "Hybrid" = {
          # Nelder-Mead
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "poi_hmm",
                           silent = TRUE)
          neldermead <- optim(par = obj$par,
                              fn = obj$fn,
                              method = "Nelder-Mead",
                              control = list(trace = 0,
                                             maxit = hybrid_iterations))
          neldermead_estimates <- pois.HMM.pw2pn(m, neldermead$par)
          working_neldermead_estimates <- pois.HMM.pn2pw(m = m,
                                                         lambda = neldermead_estimates$lambda,
                                                         gamma = neldermead_estimates$gamma)
          # nlminb on Nelder-Mead parameters
          obj <- MakeADFun(TMB_data,
                           working_neldermead_estimates,
                           DLL = "poi_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        # setup = {
        #   obj <<- MakeADFun(TMB_data,
        #                     working_params,
        #                     DLL = "poi_hmm",
        #                     silent = TRUE)
        # },
        times = BENCHMARK_NELDER_MEAD)
        
        times <- temp_microbenchmark$time / 10^9
        time_Nelder_Mead <- times[temp_microbenchmark$expr == "Nelder-Mead"]
        time_nlminb <- times[temp_microbenchmark$expr == "nlminb"]
        time_Hybrid <- times[temp_microbenchmark$expr == "Hybrid"]
        
        time_iterations_table <- rbind(time_iterations_table,
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "tinn",
                                                  method = "Nelder-Mead",
                                                  time = time_Nelder_Mead),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "tinn",
                                                  method = "nlminb",
                                                  time = time_nlminb),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "tinn",
                                                  method = "Hybrid",
                                                  time = time_Hybrid))
        
      }
      p(sprintf("Tinn done. Starting Simu3"))
      
      # simu3 --------------------------
      if (TRUE) {
        set.seed(SEED)
        # Parameters and covariates simu3 --------------------------
        m <- M_LIST_SIMU3
        simu3_true_gamma <- matrix(c(0.95, 0.05,
                                     0.15, 0.85),
                                   byrow = TRUE,
                                   nrow = m,
                                   ncol = m)
        simu3_true_mu <- seq(-5,
                             5,
                             length.out = m)
        simu3_true_sigma <- seq(1, 5, length.out = m)
        simu3_true_delta <- stat.dist(simu3_true_gamma)
        
        simu3_data <- norm.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                         mod = list(m = m,
                                                                    mu = simu3_true_mu,
                                                                    sigma = simu3_true_sigma,
                                                                    gamma = simu3_true_gamma,
                                                                    delta = simu3_true_delta),
                                                         debug_message = "simu3_data")
        # simu3_states <- simu3_data$states
        simu3_data <- simu3_data$data
        
        # simu3_nll <- norm.HMM.mllk(parvect = norm.HMM.pn2pw(m = m,
        #                                                     mu = simu3_true_mu,
        #                                                     sigma = simu3_true_sigma,
        #                                                     gamma = simu3_true_gamma,
        #                                                     delta = simu3_true_delta),
        #                            x_alias = simu3_data,
        #                            m_alias = m)
        
        # Different parameters from the true ones, to begin estimation
        gamma <- matrix(0.2 / (m - 1),
                        nrow = m,
                        ncol = m)
        diag(gamma) <- 0.8
        mu <- seq(-2 * mean(simu3_data),
                  2 * mean(simu3_data),
                  length.out = m)
        sigma <- seq(0.5 * sd(simu3_data),
                     1.5 * sd(simu3_data),
                     length.out = m)
        delta <- stat.dist(gamma)
        
        # Parameters & covariates for TMB ------------------
        working_params <- norm.HMM.pn2pw(m = m,
                                         mu = mu,
                                         sigma = sigma,
                                         gamma = gamma,
                                         delta = delta)
        # working_params_init$true_lambda <- true_lambda
        # working_params_init$simu3_true_gamma <- simu3_true_gamma
        # working_params_init$simu3_true_delta <- simu3_true_delta
        
        TMB_data <- list(x = simu3_data,
                         m = m)
        
        # Benchmark different optimization methods (ensure high enough iterations to get convergence) -----------------------------
        temp_microbenchmark <- microbenchmark("Nelder-Mead" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "norm_hmm",
                           silent = TRUE)
          optim(par = obj$par,
                fn = obj$fn,
                method = "Nelder-Mead",
                control = CONTROL_ARGS$`Nelder-Mead`)
        },
        
        "nlminb" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "norm_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        
        "Hybrid" = {
          # Nelder-Mead
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "norm_hmm",
                           silent = TRUE)
          neldermead <- optim(par = obj$par,
                              fn = obj$fn,
                              method = "Nelder-Mead",
                              control = list(trace = 0,
                                             maxit = hybrid_iterations))
          neldermead_estimates <- norm.HMM.pw2pn(m, neldermead$par)
          working_neldermead_estimates <- norm.HMM.pn2pw(m = m,
                                                         mu = neldermead_estimates$mu,
                                                         sigma = neldermead_estimates$sigma,
                                                         gamma = neldermead_estimates$gamma,
                                                         delta = neldermead_estimates$delta)
          # nlminb on Nelder-Mead parameters
          obj <- MakeADFun(TMB_data,
                           working_neldermead_estimates,
                           DLL = "norm_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        # setup = {
        #   obj <<- MakeADFun(TMB_data,
        #                     working_params,
        #                     DLL = "norm_hmm",
        #                     silent = TRUE)
        # },
        times = BENCHMARK_NELDER_MEAD)
        
        times <- temp_microbenchmark$time / 10^9
        time_Nelder_Mead <- times[temp_microbenchmark$expr == "Nelder-Mead"]
        time_nlminb <- times[temp_microbenchmark$expr == "nlminb"]
        time_Hybrid <- times[temp_microbenchmark$expr == "Hybrid"]
        
        time_iterations_table <- rbind(time_iterations_table,
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "simu3",
                                                  method = "Nelder-Mead",
                                                  time = time_Nelder_Mead),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "simu3",
                                                  method = "nlminb",
                                                  time = time_nlminb),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "simu3",
                                                  method = "Hybrid",
                                                  time = time_Hybrid))
        
      }
      p(sprintf("Simu3 done. Starting SP500"))
      
      # SP500 --------------------------
      if (TRUE) {
        set.seed(SEED)
        CONTROL_ARGS$`Nelder-Mead`$maxit <- 10000
        # Parameters and covariates SP500 --------------------------
        m <- 3
        DATA <- sp500_scaled_data
        mu <- seq(-2 * mean(DATA),
                  2 * mean(DATA),
                  length.out = m)
        sigma <- seq(0.5 * sd(DATA),
                     1.5 * sd(DATA),
                     length.out = m)
        gamma <- matrix(0.2 / (m - 1),
                        nrow = m,
                        ncol = m)
        diag(gamma) <- 0.8
        delta <- stat.dist(gamma)
        
        # Parameters & covariates for TMB ------------------
        TMB_data <- list(x = DATA,
                         m = m)
        working_params <- norm.HMM.pn2pw(m = m,
                                         mu = mu,
                                         sigma = sigma,
                                         gamma = gamma,
                                         delta = delta)
        
        # Benchmark different optimization methods (ensure high enough iterations to get convergence) -----------------------------
        temp_microbenchmark <- microbenchmark("Nelder-Mead" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "norm_hmm",
                           silent = TRUE)
          optim(par = obj$par,
                fn = obj$fn,
                method = "Nelder-Mead",
                control = CONTROL_ARGS$`Nelder-Mead`)
        },
        
        "nlminb" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "norm_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        
        "Hybrid" = {
          # Nelder-Mead
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "norm_hmm",
                           silent = TRUE)
          neldermead <- optim(par = obj$par,
                              fn = obj$fn,
                              method = "Nelder-Mead",
                              control = list(trace = 0,
                                             maxit = hybrid_iterations))
          neldermead_estimates <- norm.HMM.pw2pn(m, neldermead$par)
          working_neldermead_estimates <- norm.HMM.pn2pw(m = m,
                                                         mu = neldermead_estimates$mu,
                                                         sigma = neldermead_estimates$sigma,
                                                         gamma = neldermead_estimates$gamma,
                                                         delta = neldermead_estimates$delta)
          # nlminb on Nelder-Mead parameters
          obj <- MakeADFun(TMB_data,
                           working_neldermead_estimates,
                           DLL = "norm_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        # setup = {
        #   obj <<- MakeADFun(TMB_data,
        #                     working_params,
        #                     DLL = "norm_hmm",
        #                     silent = TRUE)
        # },
        times = BENCHMARK_NELDER_MEAD)
        
        times <- temp_microbenchmark$time / 10^9
        time_Nelder_Mead <- times[temp_microbenchmark$expr == "Nelder-Mead"]
        time_nlminb <- times[temp_microbenchmark$expr == "nlminb"]
        time_Hybrid <- times[temp_microbenchmark$expr == "Hybrid"]
        
        time_iterations_table <- rbind(time_iterations_table,
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "sp500",
                                                  method = "Nelder-Mead",
                                                  time = time_Nelder_Mead),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "sp500",
                                                  method = "nlminb",
                                                  time = time_nlminb),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "sp500",
                                                  method = "Hybrid",
                                                  time = time_Hybrid))
      }
      p(sprintf("SP500 done. Starting Soap"))
      
      # soap --------------------------
      if (TRUE) {
        set.seed(SEED)
        CONTROL_ARGS$`Nelder-Mead`$maxit <- 10000
        # Parameters and covariates hosp --------------------------
        m <- 2
        DATA <- soap_data
        lambda <- seq(quantile(DATA, 0.1, na.rm = TRUE),
                      quantile(DATA, 0.9, na.rm = TRUE),
                      length.out = m)
        gamma <- matrix(0.2 / (m - 1),
                        nrow = m,
                        ncol = m)
        diag(gamma) <- 0.8
        delta <- stat.dist(gamma)
        
        # Parameters & covariates for TMB ------------------
        working_params <- pois.HMM.pn2pw(m = m,
                                         lambda = lambda,
                                         gamma = gamma,
                                         delta = delta)
        TMB_data <- list(x = DATA,
                         m = m)
        
        # Benchmark different optimization methods (ensure high enough iterations to get convergence) -----------------------------
        temp_microbenchmark <- microbenchmark("Nelder-Mead" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "poi_hmm",
                           silent = TRUE)
          optim(par = obj$par,
                fn = obj$fn,
                method = "Nelder-Mead",
                control = CONTROL_ARGS$`Nelder-Mead`)
        },
        
        "nlminb" = {
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "poi_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        
        "Hybrid" = {
          # Nelder-Mead
          obj <- MakeADFun(TMB_data,
                           working_params,
                           DLL = "poi_hmm",
                           silent = TRUE)
          neldermead <- optim(par = obj$par,
                              fn = obj$fn,
                              method = "Nelder-Mead",
                              control = list(trace = 0,
                                             maxit = hybrid_iterations))
          neldermead_estimates <- pois.HMM.pw2pn(m, neldermead$par)
          working_neldermead_estimates <- pois.HMM.pn2pw(m = m,
                                                         lambda = neldermead_estimates$lambda,
                                                         gamma = neldermead_estimates$gamma,
                                                         delta = neldermead_estimates$delta)
          # nlminb on Nelder-Mead parameters
          obj <- MakeADFun(TMB_data,
                           working_neldermead_estimates,
                           DLL = "poi_hmm",
                           silent = TRUE)
          nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 hessian = obj$he)
        },
        # setup = {
        #   obj <<- MakeADFun(TMB_data,
        #                     working_params,
        #                     DLL = "poi_hmm",
        #                     silent = TRUE)
        # },
        times = BENCHMARK_NELDER_MEAD)
        
        times <- temp_microbenchmark$time / 10^9
        time_Nelder_Mead <- times[temp_microbenchmark$expr == "Nelder-Mead"]
        time_nlminb <- times[temp_microbenchmark$expr == "nlminb"]
        time_Hybrid <- times[temp_microbenchmark$expr == "Hybrid"]
        
        time_iterations_table <- rbind(time_iterations_table,
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "soap",
                                                  method = "Nelder-Mead",
                                                  time = time_Nelder_Mead),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "soap",
                                                  method = "nlminb",
                                                  time = time_nlminb),
                                       data.frame(iterations = hybrid_iterations,
                                                  data = "soap",
                                                  method = "Hybrid",
                                                  time = time_Hybrid))
        
      }
      p(sprintf("Soap done. Starting Hosp"))
      
    }
  })
}

# Map initial conditions SP500
if (TRUE) {
  SEED <- 1
  # CONTROL_ARGS$`Nelder-Mead`$maxit <- 10000
  set.seed(SEED)
  data <- sp500_scaled_data
  # Parameters and covariates --------------------------
  m <- 2
  gamma_init <- matrix(0.1 / (m - 1),
                       nrow = m,
                       ncol = m)
  diag(gamma_init) <- 0.9
  mu_init <- seq(-2 * mean(data),
                 2 * mean(data),
                 length.out = m)
  sigma_init <- seq(0.5 * sd(data),
                    2 * sd(data),
                    length.out = m)
  delta_init <- stat.dist(gamma_init)
  
  # Parameters & covariates for TMB ------------------
  working_params_init <- norm.HMM.pn2pw(m,
                                        mu_init,
                                        sigma_init,
                                        gamma_init,
                                        delta_init)
  
  TMB_data <- list(x = data,
                   m = m)
  
  # Estimation and MLE checks ------------------------------------
  result_BFGS <- norm.TMB.estimate(TMB_data = TMB_data,
                                   working_parameters = working_params_init,
                                   std_error = TRUE,
                                   optimizer = "BFGS")
  result_BFGS_gr <- norm.TMB.estimate(TMB_data = TMB_data,
                                      working_parameters = working_params_init,
                                      gradient = TRUE,
                                      std_error = TRUE,
                                      optimizer = "BFGS")
  result_L_BFGS_B <- norm.TMB.estimate(TMB_data = TMB_data,
                                       working_parameters = working_params_init,
                                       std_error = TRUE,
                                       optimizer = "L-BFGS-B")
  result_L_BFGS_B_gr <- norm.TMB.estimate(TMB_data = TMB_data,
                                          working_parameters = working_params_init,
                                          gradient = TRUE,
                                          std_error = TRUE,
                                          optimizer = "L-BFGS-B")
  result_CG <- norm.TMB.estimate(TMB_data = TMB_data,
                                 working_parameters = working_params_init,
                                 std_error = TRUE,
                                 optimizer = "CG")
  result_CG_gr <- norm.TMB.estimate(TMB_data = TMB_data,
                                    working_parameters = working_params_init,
                                    gradient = TRUE,
                                    std_error = TRUE,
                                    optimizer = "CG")
  result_Nelder_Mead <- norm.TMB.estimate(TMB_data = TMB_data,
                                          working_parameters = working_params_init,
                                          std_error = TRUE,
                                          optimizer = "Nelder-Mead")
  result_nlm <- norm.TMB.estimate(TMB_data = TMB_data,
                                  working_parameters = working_params_init,
                                  std_error = TRUE,
                                  optimizer = "nlm")
  result_nlm_gr <- norm.TMB.estimate(TMB_data = TMB_data,
                                     working_parameters = working_params_init,
                                     gradient = TRUE,
                                     std_error = TRUE,
                                     optimizer = "nlm")
  result_nlm_he <- norm.TMB.estimate(TMB_data = TMB_data,
                                     working_parameters = working_params_init,
                                     hessian = TRUE,
                                     std_error = TRUE,
                                     optimizer = "nlm")
  result_nlm_grhe <- norm.TMB.estimate(TMB_data = TMB_data,
                                       working_parameters = working_params_init,
                                       gradient = TRUE,
                                       hessian = TRUE,
                                       std_error = TRUE,
                                       optimizer = "nlm")
  result_nlminb <- norm.TMB.estimate(TMB_data = TMB_data,
                                     working_parameters = working_params_init,
                                     std_error = TRUE,
                                     optimizer = "nlminb")
  result_nlminb_gr <- norm.TMB.estimate(TMB_data = TMB_data,
                                        working_parameters = working_params_init,
                                        gradient = TRUE,
                                        std_error = TRUE,
                                        optimizer = "nlminb")
  result_nlminb_he <- norm.TMB.estimate(TMB_data = TMB_data,
                                        working_parameters = working_params_init,
                                        hessian = TRUE,
                                        std_error = TRUE,
                                        optimizer = "nlminb")
  result_nlminb_grhe <- norm.TMB.estimate(TMB_data = TMB_data,
                                          working_parameters = working_params_init,
                                          gradient = TRUE,
                                          hessian = TRUE,
                                          std_error = TRUE,
                                          optimizer = "nlminb")
  result_hjn <- norm.TMB.estimate(TMB_data = TMB_data,
                                  working_parameters = working_params_init,
                                  gradient = TRUE,
                                  hessian = TRUE,
                                  std_error = TRUE,
                                  optimizer = "hjn")
  result_marqLevAlg <- norm.TMB.estimate(TMB_data = TMB_data,
                                         working_parameters = working_params_init,
                                         gradient = TRUE,
                                         hessian = TRUE,
                                         std_error = TRUE,
                                         optimizer = "marqLevAlg")
  result_marqLevAlg_gr <- norm.TMB.estimate(TMB_data = TMB_data,
                                            working_parameters = working_params_init,
                                            gradient = TRUE,
                                            std_error = TRUE,
                                            optimizer = "marqLevAlg")
  result_marqLevAlg_he <- norm.TMB.estimate(TMB_data = TMB_data,
                                            working_parameters = working_params_init,
                                            hessian = TRUE,
                                            std_error = TRUE,
                                            optimizer = "marqLevAlg")
  result_marqLevAlg_grhe <- norm.TMB.estimate(TMB_data = TMB_data,
                                              working_parameters = working_params_init,
                                              gradient = TRUE,
                                              hessian = TRUE,
                                              std_error = TRUE,
                                              optimizer = "marqLevAlg")
  result_newuoa <- norm.TMB.estimate(TMB_data = TMB_data,
                                     working_parameters = working_params_init,
                                     gradient = TRUE,
                                     hessian = TRUE,
                                     std_error = TRUE,
                                     optimizer = "newuoa")
  result_BBoptim <- norm.TMB.estimate(TMB_data = TMB_data,
                                      working_parameters = working_params_init,
                                      std_error = TRUE,
                                      optimizer = "BBoptim")
  
  # If one doesn't converge successfully, stop
  # if (result_BFGS$convergence == FALSE) {
  #   stop(paste("BFGS didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_BFGS_gr$convergence == FALSE) {
  #   stop(paste("BFGS with gradient didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_L_BFGS_B$convergence == FALSE) {
  #   stop(paste("L_BFGS_B didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_L_BFGS_B_gr$convergence == FALSE) {
  #   stop(paste("L_BFGS_B with gradient didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_CG$convergence == FALSE) {
  #   stop(paste("CG didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_CG_gr$convergence == FALSE) {
  #   stop(paste("CG with gradient didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_Nelder_Mead$convergence == FALSE) {
  #   stop(paste("Nelder_Mead didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlm$convergence == FALSE) {
  #   stop(paste("nlm didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlm_gr$convergence == FALSE) {
  #   stop(paste("nlm with gradient didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlm_he$convergence == FALSE) {
  #   stop(paste("nlm with hessian didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlm_grhe$convergence == FALSE) {
  #   stop(paste("nlm with gradient + hessian didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlminb$convergence == FALSE) {
  #   stop(paste("nlminb didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlminb_gr$convergence == FALSE) {
  #   stop(paste("nlminb with gradient didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlminb_he$convergence == FALSE) {
  #   stop(paste("nlminb with hessian didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_nlminb_grhe$convergence == FALSE) {
  #   stop(paste("nlminb with gradient + hessian didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_hjn$convergence == FALSE) {
  #   stop(paste("hjn didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_marqLevAlg$convergence == FALSE) {
  #   stop(paste("marqLevAlg didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_marqLevAlg_gr$convergence == FALSE) {
  #   stop(paste("marqLevAlg with gradient didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_marqLevAlg_he$convergence == FALSE) {
  #   stop(paste("marqLevAlg with hessian didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_marqLevAlg_grhe$convergence == FALSE) {
  #   stop(paste("marqLevAlg with gradient + hessian didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_newuoa$convergence == FALSE) {
  #   stop(paste("newuoa didn't converge properly, sp500 dataset, m =", m))
  # }
  # if (result_BBoptim$convergence == FALSE) {
  #   stop(paste("BBoptim didn't converge properly, sp500 dataset, m =", m))
  # }
  
  nll_sp500 <- c("BFGS" = result_BFGS$nll,
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
                 "newuoa" = result_newuoa$nll,
                 "BBoptim" = result_BBoptim$nll)
  
  optimizers_estimates_sp500 <- list("BFGS" = result_BFGS[PARAMS_NAMES_NORM],
                                     "BFGS_gr" = result_BFGS_gr[PARAMS_NAMES_NORM],
                                     "L_BFGS_B" = result_L_BFGS_B[PARAMS_NAMES_NORM],
                                     "L_BFGS_B_gr" = result_L_BFGS_B_gr[PARAMS_NAMES_NORM],
                                     "CG" = result_CG[PARAMS_NAMES_NORM],
                                     "CG_gr" = result_CG_gr[PARAMS_NAMES_NORM],
                                     "Nelder_Mead" = result_Nelder_Mead[PARAMS_NAMES_NORM],
                                     "nlm" = result_nlm[PARAMS_NAMES_NORM],
                                     "nlm_gr" = result_nlm_gr[PARAMS_NAMES_NORM],
                                     "nlm_he" = result_nlm_he[PARAMS_NAMES_NORM],
                                     "nlm_grhe" = result_nlm_grhe[PARAMS_NAMES_NORM],
                                     "nlminb" = result_nlminb[PARAMS_NAMES_NORM],
                                     "nlminb_gr" = result_nlminb_gr[PARAMS_NAMES_NORM],
                                     "nlminb_he" = result_nlminb_he[PARAMS_NAMES_NORM],
                                     "nlminb_grhe" = result_nlminb_grhe[PARAMS_NAMES_NORM],
                                     "hjn" = result_hjn[PARAMS_NAMES_NORM],
                                     "marqLevAlg" = result_marqLevAlg[PARAMS_NAMES_NORM],
                                     "marqLevAlg_gr" = result_marqLevAlg_gr[PARAMS_NAMES_NORM],
                                     "marqLevAlg_he" = result_marqLevAlg_he[PARAMS_NAMES_NORM],
                                     "marqLevAlg_grhe" = result_marqLevAlg_grhe[PARAMS_NAMES_NORM],
                                     "newuoa" = result_newuoa[PARAMS_NAMES_NORM],
                                     "BBoptim" = result_BBoptim[PARAMS_NAMES_NORM])
  
  
  # Creating useful variables -----------------
  mu_indices <- 1:m
  sigma_indices <- max(mu_indices) + 1:m
  gamma_indices <- max(sigma_indices) + 1:(m ^ 2)
  mu_names <- paste0("mu", 1:m)
  sigma_names <- paste0("sigma", 1:m)
  gamma_names <- paste0("gamma", rep(1:m, each = m), 1:m)
  # delta_indices <- max(sigma_indices) + (1:m)
  tgamma_indices <- max(sigma_indices) + 1:(m ^ 2 - m)
  
  # mus_rowwise <- sapply(optimizers_estimates_sp500, "[[", "mu")
  # sigmas_rowwise <- sapply(optimizers_estimates_sp500, "[[", "sigma")
  # gammas_rowwise <- sapply(optimizers_estimates_sp500, "[[", "gamma")
  # deltas_rowwise <- sapply(optimizers_estimates_sp500, "[[", "delta")
  # true_mu <- apply(mus_rowwise, 1, quantile, probs = 0.75)
  # true_sigma <- apply(sigmas_rowwise, 1, quantile, probs = 0.75)
  # true_gamma <- apply(gammas_rowwise, 1, quantile, probs = 0.75)
  # true_gamma <- matrix(true_gamma, nrow = m, ncol = m)
  # true_delta <- apply(deltas_rowwise, 1, quantile, probs = 0.75)
  
  ## Create parameters variations -------------------
  sp500_sd <- sd(data)
  accuracy_mus <- seq(-20, 20, length = 10)
  accuracy_sigmas <- seq(0.1, 5 * sp500_sd, length = 10)
  accuracy_gammas <- (1:9) / 10
  
  params_mu <- setNames(lapply(1:m, function(x) bquote(accuracy_mus)),
                        mu_names)
  params_sigma <- setNames(lapply(1:m, function(x) bquote(accuracy_sigmas)),
                           sigma_names)
  params_gamma <- setNames(lapply(1:(m^2), function(x) bquote(accuracy_gammas)),
                           gamma_names)
  
  ## Generate possible mus + sigmas -------------------
  parameter_possibilities_mu <- do.call("CJ", params_mu)
  parameter_possibilities_sigma <- do.call("CJ", params_sigma)
  
  parameter_possibilities_mu$idx <- 1:nrow(parameter_possibilities_mu)
  
  # Remove rows where the mus are not strictly increasing or unique
  idx_rows_to_remove <- apply(X = parameter_possibilities_mu,
                              MARGIN = 1,
                              FUN = function(row, m2 = m) {
                                
                                row_idx <- row[["idx"]]
                                
                                row_values_mu <- row[names(row) != "idx"]
                                # Are any mus not sorted (not increasing), and are any mus not unique (not strict order)
                                if (any(duplicated(row_values_mu)) == TRUE || any(row_values_mu != sort(row_values_mu))) {
                                  return(row_idx)
                                }
                                
                              })
  # Remove NULL elements and convert to vector at the same time (vectors cannot contain NULL)
  idx_rows_to_remove <- unlist(idx_rows_to_remove)
  parameter_possibilities_mu <- parameter_possibilities_mu[-idx_rows_to_remove, ]
  
  ## Generate possible gammas ---------------------
  parameter_possibilities_gamma <- do.call("CJ", params_gamma)
  
  parameter_possibilities_gamma$idx <- 1:nrow(parameter_possibilities_gamma)
  
  # Remove rows where the TPM rows do not sum to 1
  idx_rows_to_remove <- apply(X = parameter_possibilities_gamma,
                              MARGIN = 1,
                              FUN = function(row, m2 = m) {
                                
                                row_idx <- row[["idx"]]
                                
                                row_values_gamma <- row[names(row) != "idx"]
                                # We go through each row 1:m (i.e. column indices 1:m, (m+1):(m+m), ..., (m^2-m+1):(m^2))
                                incremental_indices_gamma <- 1:m2
                                for (counter in 1:m2) {
                                  if (sum(row_values_gamma[incremental_indices_gamma]) != 1) {
                                    return(row_idx)
                                  }
                                  incremental_indices_gamma <- incremental_indices_gamma + m
                                }
                                
                              })
  # Remove NULL elements and convert to vector at the same time (vectors cannot contain NULL)
  idx_rows_to_remove <- unlist(idx_rows_to_remove)
  
  parameter_possibilities_gamma <- parameter_possibilities_gamma[-idx_rows_to_remove, ]
  
  # ## Generate stopping criterion -------------------------
  # parameter_possibilities_stopping_criterion <- cbind("stopping_criterion" = 10^(-2:-5))
  
  ## Prepare possibilities --------------------------
  # data.frame objects are easier to manipulate than data.table objects
  parameter_possibilities_mu <- as.data.frame(parameter_possibilities_mu)
  parameter_possibilities_sigma <- as.data.frame(parameter_possibilities_sigma)
  parameter_possibilities_gamma <- as.data.frame(parameter_possibilities_gamma)
  # Remove idx column
  parameter_possibilities_mu <- parameter_possibilities_mu[, names(parameter_possibilities_mu) != "idx"]
  # parameter_possibilities_sigma did not get an idx column
  parameter_possibilities_gamma <- parameter_possibilities_gamma[, names(parameter_possibilities_gamma) != "idx"]
  
  full_parameter_possibilities <- merge(x = parameter_possibilities_mu, y = parameter_possibilities_sigma)
  full_parameter_possibilities <- merge(x = full_parameter_possibilities, y = parameter_possibilities_gamma)
  # full_parameter_possibilities <- merge(x = full_parameter_possibilities, y = parameter_possibilities_stopping_criterion)
  
  sample_indices <- sample(x = 1:nrow(full_parameter_possibilities), size = SAMPLE_SIZE_INITIAL_CONDITIONS)
  full_parameter_possibilities <- full_parameter_possibilities[sample_indices, ]
  
  full_parameter_possibilities$idx <- 1:nrow(full_parameter_possibilities)
  
  # Execute all parameter possibilities ------------------
  
  # registerDoSEQ()
  registerDoFuture()
  if (Sys.info()[['sysname']] == "Windows") {
    plan(multisession, workers = CORES-1)
  } else {
    plan(multicore, workers = CORES-1)
  }
  
  set.seed(SEED)
  
  oldw <- getOption("warn")
  options(warn = -1)
  # parameter_possibilities_row <- full_parameter_possibilities[75, ]
  with_progress({
    p <- progressor(along = 1:nrow(full_parameter_possibilities))
    convergence_table <- foreach(parameter_possibilities_row = iter(full_parameter_possibilities, by = "row"),
                                 .packages = c(OPTIMIZER_PACKAGES),
                                 .options.future = list(chunk.size = 20),
                                 .combine = rbind,
                                 .inorder = TRUE) %dorng% {
                                   results <- data.frame()
                                   # for (i in 1:nrow(full_parameter_possibilities)) {
                                   #   parameter_possibilities_row <- full_parameter_possibilities[i, ]
                                   
                                   dyn.load(dynlib("code/norm_hmm"))
                                   p()
                                   
                                   # Variables ---------
                                   METHODS <- c("Nelder_Mead", "Hybrid")
                                   accuracy_rates_names <- c(paste0(METHODS, "_convergence_failure"), paste0(METHODS, "_nll_found"))
                                   accuracy_rates_names <- setNames(lapply(accuracy_rates_names, function(x) NA),
                                                                    accuracy_rates_names)
                                   accuracy_rates_sp500 <- do.call(data.frame, accuracy_rates_names)
                                   accuracy_rates_sp500 <- accuracy_rates_sp500[, order(names(accuracy_rates_sp500))]
                                   accuracy_rates_sp500 <- cbind(m = m,
                                                                 accuracy_rates_sp500,
                                                                 idx = parameter_possibilities_row$idx)
                                   
                                   
                                   # Starting parameters with some variation -----------
                                   accuracy_mu <- as.numeric(parameter_possibilities_row[, mu_indices])
                                   accuracy_sigma <- as.numeric(parameter_possibilities_row[, sigma_indices])
                                   accuracy_gamma_temp <- as.numeric(parameter_possibilities_row[, gamma_indices])
                                   accuracy_gamma <- matrix(accuracy_gamma_temp, nrow = m, ncol = m, byrow = TRUE)
                                   
                                   reference_nll_accuracy <- median(nll_sp500)
                                   
                                   # reference_mu_accuracy <- true_mu
                                   # reference_sigma_accuracy <- true_sigma
                                   # reference_gamma_accuracy <- true_gamma
                                   # reference_delta_accuracy <- true_delta
                                   
                                   # Parameters & covariates for TMB ----------
                                   accuracy_working_params <- norm.HMM.pn2pw(m, accuracy_mu, accuracy_sigma, accuracy_gamma)
                                   TMB_accuracy_data <- list(x = data,
                                                             m = m)
                                   
                                   # Optimization neldermead alone / hybrid -----------------
                                   # nll_difference <- NA
                                   # check_criterion <- NA
                                   nlminb_alone_conv <- -1
                                   Hybrid_conv <- -1
                                   suppressMessages({
                                     
                                     ## nlminb alone
                                     result_nlminb_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                                 working_parameters = accuracy_working_params,
                                                                                 gradient = TRUE,
                                                                                 hessian = TRUE,
                                                                                 optimizer = "nlminb")
                                     nlminb_alone_conv <- if (exists("result_nlminb_accuracy") && !is.null(result_nlminb_accuracy)) result_nlminb_accuracy$convergence else NA
                                     
                                     ## Hybrid: Nelder-Mead + nlminb
                                     for(iterations in c(1, seq(10, 1e4, by = 10))) {
                                       
                                       result_Nelder_Mead_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                                        working_parameters = accuracy_working_params,
                                                                                        control_list = list("Nelder-Mead" = list(trace = 0,
                                                                                                                                 maxit = iterations)),
                                                                                        optimizer = "Nelder-Mead")
                                       
                                       if (exists("result_Nelder_Mead_accuracy") && !is.null(result_Nelder_Mead_accuracy)) {
                                         result_Hybrid_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                                     working_parameters = norm.HMM.pn2pw(m = result_Nelder_Mead_accuracy$m,
                                                                                                                         mu = result_Nelder_Mead_accuracy$mu,
                                                                                                                         sigma = result_Nelder_Mead_accuracy$sigma,
                                                                                                                         gamma = result_Nelder_Mead_accuracy$gamma),
                                                                                     gradient = TRUE,
                                                                                     hessian = TRUE,
                                                                                     optimizer = "nlminb")
                                         Hybrid_conv <- if (exists("result_Hybrid_accuracy") && !is.null(result_Hybrid_accuracy)) result_Hybrid_accuracy$convergence else NA
                                       } else {
                                         Hybrid_conv <- NA
                                       }
                                       
                                       ## if Hybrid converges, end the loop
                                       if (exists("result_Hybrid_accuracy") && !is.null(result_Hybrid_accuracy) && result_Hybrid_accuracy$convergence == TRUE) {
                                         # nll_difference <- result_Hybrid_accuracy$nll - reference_nll_accuracy
                                         # check_criterion <- abs(nll_difference) < parameter_possibilities_row$stopping_criterion
                                         # if (check_criterion) break
                                         break
                                       }
                                     }
                                     
                                     ## Record results
                                     results <- rbind(results, data.frame(parameter_possibilities_row,
                                                                          iterations,
                                                                          Hybrid_conv,
                                                                          # nll_difference,
                                                                          # check_criterion,
                                                                          nlminb_alone_conv))
                                   })
                                   
                                   
                                   # # Return ------------------
                                   # results <- rbind(results, data.frame(parameter_possibilities_row,
                                   #                                      iterations,
                                   #                                      Hybrid_conv,
                                   #                                      # nll_difference,
                                   #                                      # check_criterion,
                                   #                                      nlminb_alone_conv))
                                   return(results)
                                 }
  })
  
  options(warn = oldw)
  
  # convergence_table %>%
  #   filter(check == TRUE) %>%
  #   mutate(stopping_criterion = factor(stopping_criterion)) %>%
  # ggplot(aes(x = stopping_criterion)) +
  #   geom_boxplot(aes(y = iterations)) +
  #   geom_text(stat = 'count', aes(label = after_stat(count)), position = position_stack(vjust=0.2)) +
  #   ggtitle("Iterations vs stopping criterion") +
  #   scale_x_discrete(limits = rev)
}

# save(time_iterations_table, file = "data/time_iterations_table.RData")
save(convergence_table, file = "data/neldermead_convergence_table.RData")

# beep()

# load("data/better_neldermead_convergence_table.RData")
# convergence_table %>%
#   select(c("idx", "stopping_criterion", "iterations", "Hybrid_conv", "nll_difference", "check_criterion", "nlminb_alone_conv")) %>%
#   rename(nlminb_conv = nlminb_alone_conv) %>%
#   pivot_longer(cols = c("Hybrid_conv", "nlminb_conv"),
#                names_pattern = "(.*)_conv",
#                values_to = "convergence") %>%
#   rename(algorithm = name) %>%
#   # Replace estimation failure (NA) with FALSE
#   mutate(convergence = case_when(is.na(convergence) ~ FALSE,
#                                  TRUE ~ convergence)) %>%
#   ggplot(aes(x = algorithm,
#              fill = convergence)) +
#   geom_bar() +
#   facet_grid(. ~ iterations, space ="free_x", scales="free_x", switch="x") +
#   geom_text(stat = 'count', aes(label = after_stat(count)), position = position_stack(vjust=0.5)) +
#   theme_Publication() +
#   theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) 

# 
# 
# load("data/time_iterations_table.RData")
# for (data_set in unique(time_iterations_table$data)) {
# 
#   ## BOXPLOT
# 
#   # print(time_iterations_table %>%
#   #         mutate(method = factor(method, levels = unique(method))) %>%
#   #         filter(data == data_set) %>%
#   #         ggplot(aes(x = method, y = time)) +
#   #         # geom_boxplot(aes(fill = method), show.legend = FALSE) +
#   #         geom_bar(aes(fill = method), show.legend = FALSE, stat = "summary", fun = "mean") +
#   #         facet_grid(. ~ iterations, space ="free_x", scales="free_x", switch="x") +
#   #         ggtitle(paste("Time vs iterations, ", data_set)))
# 
#   ## LINE PLOT (mean only)
#   print(time_iterations_table %>%
#           mutate(method = factor(method, levels = unique(method))) %>%
#           filter(data == data_set) %>%
#           group_by(iterations, method) %>%
#           summarise(mean_time = mean(time),
#                     min_time = min(time),
#                     max_time = max(time),
#                     .groups = "keep") %>%
#           ggplot(aes(x = iterations, y = mean_time)) +
#           geom_line(aes(colour = method)) +
#           # facet_grid(. ~ iterations, space ="free_x", scales="free_x", switch="x") +
#           ggtitle(paste("Time vs iterations, ", data_set)) +
#           ylab("Time") +
#           xlab("Iterations") +
#           theme_Publication())
# }
