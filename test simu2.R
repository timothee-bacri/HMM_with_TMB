load("data/results_simu2_1.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_2.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_3.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_4.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_5.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_6.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_7.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_8.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
load("data/results_simu2_9.RData")
print.data.frame(conf_int_simu2[1:12, c(2,3,4,11,12)], row.names = FALSE)
# prmatrix(matrix(conf_int_simu2$Coverage.TMB[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))
# prmatrix(matrix(conf_int_simu2$Coverage.Profile[4:12], nrow = 3, ncol = 3, byrow = TRUE), rowlab = rep("",3), collab = rep("",3))



# setwd(dir = "..")
source("code/setup_parameters.R") # NEED TO COMMENT rm(list=ls()) IN ORDER TO RUN AS A LOCAL JOB

true_lambda <- seq(1, 7, length.out = 3)
# true_lambda <- seq(1, 11, length.out = 3)
# true_lambda <- seq(1, 20, length.out = 3)

true_gamma <- matrix(c(0.95, 0.025, 0.025,
                       0.05, 0.9, 0.05,
                       0.075, 0.075, 0.85), byrow = TRUE, nrow = 3, ncol = 3)

# true_gamma <- matrix(c(0.9, 0.05, 0.05,
#                        0.05, 0.9, 0.05,
#                        0.075, 0.075, 0.85), byrow = TRUE, nrow = 3, ncol = 3)

# true_gamma <- matrix(c(0.95, 0.025, 0.025,
#                        0.05, 0.9, 0.05,
#                        0.05, 0.05, 0.9), byrow = TRUE, nrow = 3, ncol = 3)
begin <- Sys.time()
idx=1
set.seed(123)
# Parameters and covariates --------------------------
m <- M_LIST_SIMU2[idx]
# if (m == 1) {
#   true_gamma <- matrix(1)
# } else {
#   true_gamma <- matrix(0.2 / (m - 1), nrow = m, ncol = m)
#   diag(true_gamma) <- 0.8
# }
# if (m == 2) {
#   true_gamma <- matrix(c(0.95, 0.05,
#                          0.15, 0.85), byrow = TRUE, nrow = m, ncol = m)
# } else if (m == 3) {
#   true_gamma <- matrix(c(0.95, 0.025, 0.025,
#                          0.05, 0.90, 0.05,
#                          0.075, 0.075, 0.85), byrow = TRUE, nrow = m, ncol = m)
# }
# true_lambda <- seq(1, 7, length.out = m)
true_delta <- stat.dist(true_gamma)

simu2_data <- pois.HMM.generate_sample(ns = DATA_SIZE_SIMU2,
                                       mod = list(m = m,
                                                  lambda = true_lambda,
                                                  gamma = true_gamma,
                                                  delta = true_delta))$data

# Different parameters from the true ones, to begin estimation
if (m != 1) {
  gamma_init <- matrix(0.1 / (m - 1), nrow = m, ncol = m)
  diag(gamma_init) <- 0.9
} else {
  gamma_init <- matrix(1)
}
lambda_init <- seq(2, 6, length.out = m)
delta_init <- stat.dist(gamma_init)

# Parameters & covariates for TMB ------------------
working_params_init <- pois.HMM.pn2pw(m, lambda_init, gamma_init)
TMB_data <- list(x = simu2_data, m = m)
obj_init <- MakeADFun(TMB_data, working_params_init, DLL = "poi_hmm", silent = TRUE)
parvect_init <- pois.HMM.pn2pw(m = m, lambda = lambda_init, gamma = gamma_init, delta = delta_init)
parvect_init <- unlist(parvect_init)

# Estimation ------------------------------------
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
if (tmb$convergence != 0) {
  stop(paste("tmb didn't converge properly, simu2 dataset, m =", m))
}
if (tmb_g$convergence != 0) {
  stop(paste("tmb_g didn't converge properly, simu2 dataset, m =", m))
}
if (tmb_h$convergence != 0) {
  stop(paste("tmb_h didn't converge properly, simu2 dataset, m =", m))
}
if (tmb_gh$convergence != 0) {
  stop(paste("tmb_gh didn't converge properly, simu2 dataset, m =", m))
}

# Creating variables for the CIs -----------------
tmb_CI <- tmb_gh
params_names_latex <- paste0(rep("$\\lambda_{", m), 1:m, "}$")
for (gamma_idx in 1:m ^ 2) {
  # Get row and column indices for gamma instead of the default
  # columnwise index: the default indices are 1:m for the 1st column,
  # then (m + 1):(2 * m) for the 2nd, etc...
  row_col_idx <- matrix.col.idx.to.rowcol(gamma_idx, m)
  params_names_latex <- c(params_names_latex,
                          paste0("$\\gamma_{", toString(row_col_idx), "}$"))
}
params_names_latex <- c(params_names_latex,
                        paste0(rep("$\\delta_{", m), 1:m, "}$"))
len_par <- length(params_names_latex)
indices <- (length(conf_int_simu2$m) + 1):(length(conf_int_simu2$m) + len_par)
conf_int_simu2[indices, "m"] <- m
conf_int_simu2[indices, "Parameter"] <- params_names_latex
# Reminder, PARAMS_NAMES contains c("lambda", "gamma", "delta")
conf_int_simu2[indices, "Estimate"] <- unlist(tmb_CI[PARAMS_NAMES])
conf_int_simu2[indices, "True.value"] <- as.numeric(c(true_lambda, true_gamma, true_delta))

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

# TMB confidence intervals --------------
# Manually cap values at their natural bound
# lambda must be strictly above 0
lambda_L <- pmax(0.0001, tmb_gh$lambda - q95_norm * tmb_gh$lambda_std_error)
# gamma must be 0 or more
gamma_L <- pmax(0, tmb_gh$gamma - q95_norm * tmb_gh$gamma_std_error)
# delta must be 0 or above
delta_L <- pmax(0, tmb_gh$delta - q95_norm * tmb_gh$delta_std_error)
conf_int_simu2$TMB.L[which(conf_int_simu2$m == m)] <- c(lambda_L,
                                                        gamma_L,
                                                        delta_L)
# no upper bound on lambda
# gamma must be 1 or less
gamma_U <- pmin(1, tmb_gh$gamma + q95_norm * tmb_gh$gamma_std_error)
# delta must be 1 or less
delta_U <- pmin(1, tmb_gh$delta + q95_norm * tmb_gh$delta_std_error)
conf_int_simu2$TMB.U[which(conf_int_simu2$m == m)] <- c(tmb_gh$lambda + q95_norm * tmb_gh$lambda_std_error,
                                                        gamma_U,
                                                        delta_U)
# Coverage probabilities of the 3 CI methods -----------------
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
  # Unlike with the other datasets, we know the true parameters of this one
  coverage_model <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_SIMU2,
                                                       mod = list(m = m,
                                                                  lambda = true_lambda,
                                                                  gamma = true_gamma,
                                                                  delta = true_delta),
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
    working_conf_int$lower[working_conf_int$w_parameter == w_param_name] <- ci[1]
    working_conf_int$upper[working_conf_int$w_parameter == w_param_name] <- ci[2]
  }
  
  lambda_profile_lower <- exp(working_conf_int$lower[lambda_indices])
  lambda_profile_upper <- exp(working_conf_int$upper[lambda_indices])
  gamma_profile_lower <- as.numeric(Gamma_w2n(m, working_conf_int$lower[tgamma_indices]))
  gamma_profile_upper <- as.numeric(Gamma_w2n(m, working_conf_int$upper[tgamma_indices]))
  
  print(lambda_profile_lower)
  print(lambda_profile_upper)
  print(gamma_profile_lower)
  print(gamma_profile_upper)
  # If profiling doesn't yield results for all parameters, try a new coverage sample
  if (anyNA(c(lambda_profile_lower, lambda_profile_upper, gamma_profile_lower, gamma_profile_upper), recursive = TRUE)) {
    idx_coverage <- idx_coverage - 1
    profile_skips_simu2[idx] <- profile_skips_simu2[idx] + 1
    print(profile_skips_simu2[idx])
    next
  }
  # If the true value for lambda is in the CI, then increase the count
  real_lambda_profile_lower <- pmin(lambda_profile_lower, lambda_profile_upper)
  real_lambda_profile_upper <- pmax(lambda_profile_lower, lambda_profile_upper)
  indices <- which(true_lambda >= real_lambda_profile_lower & true_lambda <= real_lambda_profile_upper)
  coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
  
  # Same for gamma
  real_gamma_profile_lower <- pmin(gamma_profile_lower, gamma_profile_upper)
  real_gamma_profile_upper <- pmax(gamma_profile_lower, gamma_profile_upper)
  indices <- which(as.vector(true_gamma) >= real_gamma_profile_lower & as.vector(true_gamma) <= real_gamma_profile_upper)
  indices <- indices + m
  coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
  
  # # Sometimes, there are no profile CIs for gamma because of the need to transform the CIs to their natural form
  # # In these cases, we reduce the total coverage_adjusted_sample_size_gamma_simu2[idx, "size"]
  # # Either or both bounds for gamma can be missing due to missing tgamma values
  # # If there are missing tgamma values, the confidence interval has 0 or 1 bound
  # # and is thus useless for coverage probabilities
  # if (anyNA(gamma_profile_lower, recursive = TRUE) || anyNA(gamma_profile_upper, recursive = TRUE)) {
  #   coverage_adjusted_sample_size_gamma_simu2[idx, "size"] <- coverage_adjusted_sample_size_gamma_simu2[idx, "size"] - 1
  # } else {
  #   # Same as above for gamma
  #   indices <- which(as.vector(tmb_CI$gamma) >= gamma_profile_lower & as.vector(tmb_CI$gamma) <= gamma_profile_upper)
  #   indices <- indices + m
  #   coverage_count_profile[indices, "count"] <- coverage_count_profile[indices, "count"] + 1
  # }
  
  # Confidence interval bootstrap -----------------------------------
  bootstrap_simu2 <- data.frame()
  lambda_coverage <- coverage_model$mod$lambda
  gamma_coverage <- coverage_model$mod$gamma
  delta_coverage <- coverage_model$mod$delta
  if (BOOTSTRAP_SAMPLES != 0) {
    for (idx_sample in 1:BOOTSTRAP_SAMPLES) {
      temp <- pois.HMM.generate_estimable_sample(ns = DATA_SIZE_SIMU2,
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
      bootstrap_simu2[idx_sample, 1:len_par] <- natural_parameters
    }
    # colnames(bootstrap_simu2) <- params_names_latex
    q <- apply(bootstrap_simu2, 2, quantile.colwise)
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
  
}
coverage_count_profile[lambda_indices, "ratio"] <- coverage_count_profile[lambda_indices, "count"] / COVERAGE_SAMPLES
# coverage_count_profile[gamma_indices, "ratio"] <- coverage_count_profile[gamma_indices, "count"] / coverage_adjusted_sample_size_gamma_simu2[idx, "size"]
coverage_count_profile[gamma_indices, "ratio"] <- coverage_count_profile[gamma_indices, "count"] / COVERAGE_SAMPLES
coverage_count_profile[delta_indices, "ratio"] <- NA # delta is not a parameter for us, so it has no profile CI

coverage_count_tmb$ratio <- coverage_count_tmb$count / COVERAGE_SAMPLES

coverage_count_bootstrap$ratio <- coverage_count_bootstrap$count / COVERAGE_SAMPLES

conf_int_simu2[conf_int_simu2$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Profile"] <- coverage_count_profile$ratio * 100
conf_int_simu2[conf_int_simu2$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.Bootstrap"] <- coverage_count_bootstrap$ratio * 100
conf_int_simu2[conf_int_simu2$m == m, ][1:(m ^ 2 + 2 * m), "Coverage.TMB"] <- coverage_count_tmb$ratio * 100

# Fixes -------------------------
# Fix label switching in conf_int_simu2
ordered_params <- pois.HMM.label.order(m = m,
                                       lambda = true_lambda,
                                       gamma = true_gamma,
                                       delta = true_delta)

new_lambda_indices <- ordered_params$ordered_lambda_indices
new_gamma_indices <- ordered_params$ordered_gamma_vector_indices
new_delta_indices <- ordered_params$ordered_delta_indices

conf_int_simu2[conf_int_simu2$m == m, - 2][lambda_indices, ] <- conf_int_simu2[conf_int_simu2$m == m, - 2][lambda_indices, ][new_lambda_indices, ]
conf_int_simu2[conf_int_simu2$m == m, - 2][gamma_indices, ] <- conf_int_simu2[conf_int_simu2$m == m, - 2][gamma_indices, ][new_gamma_indices, ]
conf_int_simu2[conf_int_simu2$m == m, - 2][delta_indices, ] <- conf_int_simu2[conf_int_simu2$m == m, - 2][delta_indices, ][new_delta_indices, ]

# Reorder the TPM row-wise instead of column-wise
# Lexicographical parameter sort for gamma (sort on the parameter name)
new_gamma_indices_truncated_table <- order(conf_int_simu2[conf_int_simu2$m == m, ][gamma_indices, "Parameter"])
# Replace rows by sorted rows
conf_int_simu2[conf_int_simu2$m == m, ][gamma_indices, ] <- conf_int_simu2[gamma_indices, ][new_gamma_indices_truncated_table, ]

# The profile CIs may not be sorted, so we sort them manually
for (i in 1:length(conf_int_simu2[, 1])) {
  row <- conf_int_simu2[i, c("Profile.L", "Profile.U")]
  conf_int_simu2[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
}
mllk_values_simu2$m <- as.factor(as.numeric(mllk_values_simu2$m))
mllk_values_simu2$mllk <- as.numeric(mllk_values_simu2$mllk)
mllk_values_simu2$AIC <- as.numeric(mllk_values_simu2$AIC)
mllk_values_simu2$BIC <- as.numeric(mllk_values_simu2$BIC)
conf_int_simu2$m <- as.integer(conf_int_simu2$m)

# Remove gamma and delta when m = 1 because they're pointless
if (1 %in% M_LIST_SIMU2) {
  conf_int_simu2 <- conf_int_simu2[- c(2,3), ]
}

estim_benchmarks_df_simu2$m <- factor(estim_benchmarks_df_simu2$m, levels = M_LIST_SIMU2)
# END ---------------
simu2time <- Sys.time()
save(profile_skips_simu2, conf_int_simu2,
     file = "data/results_simu2_9.RData")

sms(paste("simu2 ", round(simu2time - begin, 1), units(simu2time - begin)))
