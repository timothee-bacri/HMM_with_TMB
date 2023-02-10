unlink("full_log_simu3.txt")
sink("full_log_simu3.txt")
idx_counter_benchmark <- 1
idx_counter_accuracy <- 1
idx_counter_bootstrap <- 1
SEED <- 1
# Prepare the data & parameters, then estimate for different numbers of hidden states
set.seed(SEED)
# Parameters and covariates --------------------------
m <- M_LIST_SIMU3
if (m == 2) {
  simu3_true_gamma <- matrix(c(0.95, 0.05,
                               0.15, 0.85),
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
simu3_true_sigma <- seq(1, 5, length.out = m)
simu3_true_delta <- stat.dist(simu3_true_gamma)

simu3_data <- norm.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                 mod = list(m = m,
                                                            mu = simu3_true_mu,
                                                            sigma = simu3_true_sigma,
                                                            gamma = simu3_true_gamma,
                                                            delta = simu3_true_delta),
                                                 debug_message = "simu3_data")
simu3_states <- simu3_data$states
simu3_data <- simu3_data$data

# Different parameters from the true ones, to begin estimation
if (m != 1) {
  gamma_init <- matrix(0.1 / (m - 1),
                       nrow = m,
                       ncol = m)
  diag(gamma_init) <- 0.9
} else {
  gamma_init <- matrix(1)
}
mu_init <- seq(-5,
               5,
               length.out = m)
sigma_init <- seq(2, 5, length.out = m)
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
# obj_init <- MakeADFun(TMB_data,
#                       working_params_init,
#                       DLL = "norm_hmm",
#                       silent = TRUE)

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
# result_ucminf <- norm.TMB.estimate(TMB_data = TMB_data,
#                               working_parameters = working_params_init,
#                               gradient = TRUE,
#                               hessian = TRUE,
#                               std_error = TRUE,
#                               optimizer = "ucminf")
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
# result_BBoptim_gr <- norm.TMB.estimate(TMB_data = TMB_data,
#                                   working_parameters = working_params_init,
#                                   gradient = TRUE,
#                                   std_error = TRUE,
#                                   optimizer = "BBoptim")


# If one doesn't converge successfully, stop
if (result_BFGS$convergence == FALSE) {
  stop(paste("BFGS didn't converge properly, simu3 dataset, m =", m))
}
if (result_BFGS_gr$convergence == FALSE) {
  stop(paste("BFGS with gradient didn't converge properly, simu3 dataset, m =", m))
}
if (result_L_BFGS_B$convergence == FALSE) {
  stop(paste("L_BFGS_B didn't converge properly, simu3 dataset, m =", m))
}
if (result_L_BFGS_B_gr$convergence == FALSE) {
  stop(paste("L_BFGS_B with gradient didn't converge properly, simu3 dataset, m =", m))
}
if (result_CG$convergence == FALSE) {
  stop(paste("CG didn't converge properly, simu3 dataset, m =", m))
}
if (result_CG_gr$convergence == FALSE) {
  stop(paste("CG with gradient didn't converge properly, simu3 dataset, m =", m))
}
if (result_Nelder_Mead$convergence == FALSE) {
  stop(paste("Nelder_Mead didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlm$convergence == FALSE) {
  stop(paste("nlm didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlm_gr$convergence == FALSE) {
  stop(paste("nlm with gradient didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlm_he$convergence == FALSE) {
  stop(paste("nlm with hessian didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlm_grhe$convergence == FALSE) {
  stop(paste("nlm with gradient + hessian didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlminb$convergence == FALSE) {
  stop(paste("nlminb didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlminb_gr$convergence == FALSE) {
  stop(paste("nlminb with gradient didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlminb_he$convergence == FALSE) {
  stop(paste("nlminb with hessian didn't converge properly, simu3 dataset, m =", m))
}
if (result_nlminb_grhe$convergence == FALSE) {
  stop(paste("nlminb with gradient + hessian didn't converge properly, simu3 dataset, m =", m))
}
if (result_hjn$convergence == FALSE) {
  stop(paste("hjn didn't converge properly, simu3 dataset, m =", m))
}
if (result_marqLevAlg$convergence == FALSE) {
  stop(paste("marqLevAlg didn't converge properly, simu3 dataset, m =", m))
}
if (result_marqLevAlg_gr$convergence == FALSE) {
  stop(paste("marqLevAlg with gradient didn't converge properly, simu3 dataset, m =", m))
}
if (result_marqLevAlg_he$convergence == FALSE) {
  stop(paste("marqLevAlg with hessian didn't converge properly, simu3 dataset, m =", m))
}
if (result_marqLevAlg_grhe$convergence == FALSE) {
  stop(paste("marqLevAlg with gradient + hessian didn't converge properly, simu3 dataset, m =", m))
}
# if (result_ucminf$convergence == FALSE) {
#   stop(paste("ucminf didn't converge properly, simu3 dataset, m =", m))
# }
if (result_newuoa$convergence == FALSE) {
  stop(paste("newuoa didn't converge properly, simu3 dataset, m =", m))
}
if (result_BBoptim$convergence == FALSE) {
  stop(paste("BBoptim didn't converge properly, simu3 dataset, m =", m))
}
# if (result_BBoptim_gr$convergence == FALSE) {
#   stop(paste("BBoptim with gradient didn't converge properly, simu3 dataset, m =", m))
# }

nll_simu3 <- c("BFGS" = result_BFGS$nll,
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


# optimizers_estimates_simu3 <- do.call("list", setNames(lapply(paste0("par_", OPTIMIZERS_METHOD_UNDERSCORE),
#                                                              function(x) bquote(pois.HMM.label.order(m = m,
#                                                                                                      lambda = get(.(x))$lambda,
#                                                                                                      gamma = get(.(x))$gamma,
#                                                                                                      delta = get(.(x))$delta,
#                                                                                                      indices = FALSE))),
#                                                       OPTIMIZERS_METHOD))
optimizers_estimates_simu3 <- list("BFGS" = result_BFGS[PARAMS_NAMES_NORM],
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
                                   # "ucminf" = result_ucminf[PARAMS_NAMES_NORM],
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

# If the estimates have a coefficient of variation too high (we pick 20%), stop
mus_rowwise <- sapply(optimizers_estimates_simu3, "[[", "mu")
sigmas_rowwise <- sapply(optimizers_estimates_simu3, "[[", "sigma")
gammas_rowwise <- sapply(optimizers_estimates_simu3, "[[", "gamma")
deltas_rowwise <- sapply(optimizers_estimates_simu3, "[[", "delta")
cv_mus <- apply(mus_rowwise, 1, sd) / rowMeans(mus_rowwise)
cv_sigmas <- apply(sigmas_rowwise, 1, sd) / rowMeans(sigmas_rowwise)
cv_gammas <- apply(gammas_rowwise, 1, sd) / rowMeans(gammas_rowwise)
cv_nll <- sd(nll_simu3) / mean(nll_simu3)
if (any(cv_mus > 0.2) || any(cv_sigmas > 0.2) || any(cv_gammas > 0.2) || any(cv_nll > 0.2)) {
  stop("Simu3 data MLEs or nlls have a coefficient of variation above 20%")
}

# Benchmarks ------------
if (BENCHMARK_SAMPLES != 0) {
  begin_benchmark <- Sys.time()
  for (idx_counter_benchmark in 1:BENCHMARK_SAMPLES) {
    set.seed(idx_counter_benchmark)
    # Generate data that can be estimated by TMB_GH
    # and is tested on the slightly off parameters from the beginning of this file
    # The goal is to have a dataset that poses no estimation problem,
    # when estimated with guessed initial parameters
    benchmark_model <- norm.HMM.generate.estimable.sample(ns = DATA_SIZE_SIMU3,
                                                          mod = list(m = m,
                                                                     mu = simu3_true_mu,
                                                                     sigma = simu3_true_sigma,
                                                                     gamma = simu3_true_gamma,
                                                                     delta = simu3_true_delta),
                                                          debug_message = "Simu3 benchmarking")
    benchmark_data <- benchmark_model$data
    # Parameters & covariates for DM and TMB
    TMB_benchmark_data <- list(x = benchmark_data,
                               m = m)
    working_params <- norm.HMM.pn2pw(m = m,
                                     mu = simu3_true_mu,
                                     sigma = simu3_true_sigma,
                                     gamma = simu3_true_gamma,
                                     delta = simu3_true_delta)
    obj_benchmark <- MakeADFun(TMB_benchmark_data,
                               working_params,
                               DLL = "norm_hmm",
                               silent = TRUE)
    
    # Benchmark different optimization methods ----------------------------------------------
    temp_microbenchmark <- microbenchmark("BFGS" = optim(par = obj_benchmark$par,
                                                         fn = obj_benchmark$fn,
                                                         method = "BFGS",
                                                         control = CONTROL_ARGS$BFGS)$convergence==0,
                                          "BFGS_gr" = optim(par = obj_benchmark$par,
                                                            fn = obj_benchmark$fn,
                                                            gr = obj_benchmark$gr,
                                                            method = "BFGS",
                                                            control = CONTROL_ARGS$BFGS)$convergence==0,
                                          "L-BFGS-B" = optim(par = obj_benchmark$par,
                                                             fn = obj_benchmark$fn,
                                                             method = "L-BFGS-B",
                                                             control = CONTROL_ARGS$`L-BFGS-B`)$convergence==0,
                                          "L-BFGS-B_gr" = optim(par = obj_benchmark$par,
                                                                fn = obj_benchmark$fn,
                                                                gr = obj_benchmark$gr,
                                                                method = "L-BFGS-B",
                                                                control = CONTROL_ARGS$`L-BFGS-B`)$convergence==0,
                                          "CG" = optim(par = obj_benchmark$par,
                                                       fn = obj_benchmark$fn,
                                                       method = "CG",
                                                       control = CONTROL_ARGS$CG)$convergence==0,
                                          "CG_gr" = optim(par = obj_benchmark$par,
                                                          fn = obj_benchmark$fn,
                                                          gr = obj_benchmark$gr,
                                                          method = "CG",
                                                          control = CONTROL_ARGS$CG)$convergence==0,
                                          "Nelder-Mead" = optim(par = obj_benchmark$par,
                                                                fn = obj_benchmark$fn,
                                                                method = "Nelder-Mead",
                                                                control = CONTROL_ARGS$`Nelder-Mead`)$convergence==0,
                                          "nlm" = nlm(f = nlm_gradient_hessian_objective,
                                                      p = obj_benchmark$par,
                                                      fn = obj_benchmark$fn,
                                                      iterlim = CONTROL_ARGS$nlm$iterlim)$code %in% 1:2,
                                          "nlm_gr" = nlm(f = nlm_gradient_hessian_objective,
                                                         p = obj_benchmark$par,
                                                         fn = obj_benchmark$fn,
                                                         gr = obj_benchmark$gr,
                                                         iterlim = CONTROL_ARGS$nlm$iterlim)$code %in% 1:2,
                                          "nlm_he" = nlm(f = nlm_gradient_hessian_objective,
                                                         p = obj_benchmark$par,
                                                         fn = obj_benchmark$fn,
                                                         he = obj_benchmark$he,
                                                         iterlim = CONTROL_ARGS$nlm$iterlim)$code %in% 1:2,
                                          "nlm_grhe" = nlm(f = nlm_gradient_hessian_objective,
                                                           p = obj_benchmark$par,
                                                           fn = obj_benchmark$fn,
                                                           gr = obj_benchmark$gr,
                                                           he = obj_benchmark$he,
                                                           iterlim = CONTROL_ARGS$nlm$iterlim)$code %in% 1:2,
                                          "nlminb" = nlminb(start = obj_benchmark$par,
                                                            objective = obj_benchmark$fn,
                                                            control = CONTROL_ARGS$nlminb)$convergence==0,
                                          "nlminb_gr" = nlminb(start = obj_benchmark$par,
                                                               objective = obj_benchmark$fn,
                                                               gradient = obj_benchmark$gr,
                                                               control = CONTROL_ARGS$nlminb)$convergence==0,
                                          "nlminb_he" = nlminb(start = obj_benchmark$par,
                                                               objective = obj_benchmark$fn,
                                                               hessian = obj_benchmark$he,
                                                               control = CONTROL_ARGS$nlminb)$convergence==0,
                                          "nlminb_grhe" = nlminb(start = obj_benchmark$par,
                                                                 objective = obj_benchmark$fn,
                                                                 gradient = obj_benchmark$gr,
                                                                 hessian = obj_benchmark$he,
                                                                 control = CONTROL_ARGS$nlminb)$convergence==0,
                                          "hjn" = hjn(par = obj_benchmark$par,
                                                      fn = obj_benchmark$fn)$convergence==0,
                                          "marqLevAlg" = marqLevAlg(b = obj_benchmark$par,
                                                                    fn = obj_benchmark$fn,
                                                                    maxiter = CONTROL_ARGS$marqLevAlg$maxiter)$istop==1,
                                          "marqLevAlg_gr" = marqLevAlg(b = obj_benchmark$par,
                                                                       fn = obj_benchmark$fn,
                                                                       gr = obj_benchmark$gr,
                                                                       maxiter = CONTROL_ARGS$marqLevAlg$maxiter)$istop==1,
                                          "marqLevAlg_he" = marqLevAlg(b = obj_benchmark$par,
                                                                       fn = obj_benchmark$fn,
                                                                       hess = obj_benchmark$he,
                                                                       maxiter = CONTROL_ARGS$marqLevAlg$maxiter)$istop==1,
                                          "marqLevAlg_grhe" = marqLevAlg(b = obj_benchmark$par,
                                                                         fn = obj_benchmark$fn,
                                                                         gr = obj_benchmark$gr,
                                                                         hess = obj_benchmark$he,
                                                                         maxiter = CONTROL_ARGS$marqLevAlg$maxiter)$istop==1,
                                          "newuoa" = newuoa(par = obj_benchmark$par,
                                                            fn = obj_benchmark$fn)$ierr==0,
                                          # "ucminf" = ucminf(par = obj_benchmark$par,
                                          #                   fn = obj_benchmark$fn,
                                          #                   gr = obj_benchmark$gr,
                                          #                   control = CONTROL_ARGS$ucminf)$convergence %in% c(1, 2, 4),
                                          "BBoptim" = BBoptim(par = obj_benchmark$par,
                                                              fn = obj_benchmark$fn,
                                                              quiet = TRUE,
                                                              control = CONTROL_ARGS$BBoptim)$convergence==0,
                                          # "BBoptim_gr" = BBoptim(par = obj_benchmark$par,
                                          #                        fn = obj_benchmark$fn,
                                          #                        gr = obj_benchmark$gr,
                                          #                        quiet = TRUE,
                                          #                        control = CONTROL_ARGS$BBoptim)$convergence==0,
                                          times = 1,
                                          check = "equal",
                                          setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data,
                                                                              working_params,
                                                                              DLL = "norm_hmm",
                                                                              silent = TRUE))
    
    # MICROBENCHMARK NLOPTR ALGORITHMS
    # nloptr_results <- data.frame(algo = character(),
    #                              status = integer(),
    #                              objective = numeric(),
    #                              message = character(),
    #                              xtol_rel = numeric(),
    #                              tlambda_min = numeric(),
    #                              iterations = integer())
    # 
    # idx_nloptr <- 1
    # for (algo in NLOPTR_ALGORITHMS) {
    #   # NEED local_opts FOR THESE ALGORITHMS
    #   if (algo %in% c("NLOPT_GN_MLSL",
    #                   "NLOPT_GD_MLSL",
    #                   "NLOPT_GN_MLSL_LDS",
    #                   "NLOPT_GD_MLSL_LDS",
    #                   "NLOPT_LN_AUGLAG",
    #                   "NLOPT_LD_AUGLAG",
    #                   "NLOPT_LN_AUGLAG_EQ",
    #                   "NLOPT_LD_AUGLAG_EQ")) {
    #     local_opts <- list(algorithm = "NLOPT_LD_MMA",
    #                        xtol_rel = 1e-10)
    #     # SOME ALGORITHMS SUCCEED WITH DIFFERENT xtol_rel. Info with nloptr.print.options():
    #     for (xtol_rel in 10^(-3:-20)) {
    #       temp_nloptr <- nloptr(x0 = obj_benchmark$par,
    #                      eval_f = obj_benchmark$fn,
    #                      ... = NULL,
    #                      eval_grad_f = obj_benchmark$gr,
    #                      opts = list(algorithm = algo,
    #                                  maxeval = 10000,
    #                                  maxtime = 120,
    #                                  xtol_rel = xtol_rel,
    #                                  ranseed = NLOPTR_SEED,
    #                                  local_opts = local_opts))
    #       # STOP WHEN THE ALGORITHM WORKS
    #       if (temp_nloptr$status == 1) {
    #         break
    #       }
    #     }
    #   } else {
    #     # SOME ALGORITHMS SUCCEED WITH DIFFERENT xtol_rel. Info with nloptr.print.options():
    #     for (xtol_rel in 10^(-3:-20)) {
    #       temp_nloptr <- nloptr(x0 = obj_benchmark$par,
    #                      eval_f = obj_benchmark$fn,
    #                      ... = NULL,
    #                      eval_grad_f = obj_benchmark$gr,
    #                      opts = list(algorithm = algo,
    #                                  maxeval = 10000,
    #                                  maxtime = 120,
    #                                  ranseed = NLOPTR_SEED,
    #                                  xtol_rel = xtol_rel))
    #       # STOP WHEN THE ALGORITHM WORKS
    #       if (temp_nloptr$status == 1) {
    #         break
    #       }
    #     }
    #   }
    #   
    #   if (temp_nloptr$status == 1) {
    #     nloptr_results[idx_nloptr, "algo"] <- algo
    #     nloptr_results[idx_nloptr, "status"] <- temp_nloptr$status
    #     nloptr_results[idx_nloptr, "objective"] <- temp_nloptr$objective
    #     nloptr_results[idx_nloptr, "message"] <- temp_nloptr$message
    #     nloptr_results[idx_nloptr, "xtol_rel"] <- xtol_rel
    #     # A list would be better to return all parameters, but it is not needed so we keep things simple.
    #     # We want to compare the lowest lambda_1 later.
    #     # Since our transformation function for the lambdas is increasing, the min tlambda corresponds to the min lambda.
    #     nloptr_results[idx_nloptr, "tlambda_min"] <- min(temp_nloptr$solution[1:m])
    #     nloptr_results[idx_nloptr, "iterations"] <- temp_nloptr$iterations
    #     
    #     idx_nloptr <- idx_nloptr + 1
    #   }
    # }
    # 
    # algos_and_xtol_rels_list <- split(nloptr_results, seq(nrow(nloptr_results)))
    # 
    # microbenchmark_params <- setNames(
    #   object = lapply(X = algos_and_xtol_rels_list,
    #                   obj_benchmark = obj_benchmark,
    #                   local_opts = local_opts,
    #                   NLOPTR_SEED = NLOPTR_SEED,
    #                   FUN = function(x, obj_benchmark, local_opts, NLOPTR_SEED) {
    #                     
    #                     if (x$algo %in% c("NLOPT_GN_MLSL",
    #                                       "NLOPT_GD_MLSL",
    #                                       "NLOPT_GN_MLSL_LDS",
    #                                       "NLOPT_GD_MLSL_LDS",
    #                                       "NLOPT_LN_AUGLAG",
    #                                       "NLOPT_LD_AUGLAG",
    #                                       "NLOPT_LN_AUGLAG_EQ",
    #                                       "NLOPT_LD_AUGLAG_EQ")) {
    #                       bquote(nloptr(x0 = .(obj_benchmark$par),
    #                                     eval_f = .(obj_benchmark$fn),
    #                                     ... = NULL,
    #                                     eval_grad_f = .(obj_benchmark$gr),
    #                                     opts = list(algorithm = .(x$algo),
    #                                                 maxeval = 10000,
    #                                                 maxtime = 120,
    #                                                 xtol_rel = .(x$xtol_rel),
    #                                                 ranseed = .(NLOPTR_SEED),
    #                                                 local_opts = .(local_opts)))$status==1)
    #                       
    #                     } else {
    #                       bquote(nloptr(x0 = .(obj_benchmark$par),
    #                                     eval_f = .(obj_benchmark$fn),
    #                                     ... = NULL,
    #                                     eval_grad_f = .(obj_benchmark$gr),
    #                                     opts = list(algorithm = .(x$algo),
    #                                                 maxeval = 10000,
    #                                                 maxtime = 120,
    #                                                 ranseed = .(NLOPTR_SEED),
    #                                                 xtol_rel = .(x$xtol_rel)))$status==1)
    #                     }
    #                   }
    #   ),
    #   nm = nloptr_results$algo
    # )
    # 
    # temp_microbenchmark_nloptr <- do.call(what = "microbenchmark",
    #                                       args = list(list = microbenchmark_params,
    #                                                   times = 1,
    #                                                   setup = obj_benchmark <<- MakeADFun(TMB_benchmark_data,
    #                                                                                       working_params_init,
    #                                                                                       DLL = "norm_hmm",
    #                                                                                       silent = TRUE)))
    
    times <- temp_microbenchmark$time / 10^9
    
    # time_list <- do.call(what = "list",
    #                      args = setNames(lapply(OPTIMIZERS_METHOD,
    #                                             function(x) bquote(times[temp_microbenchmark$expr == .(x)])),
    #                                      OPTIMIZERS_METHOD))
    time_BFGS <- times[temp_microbenchmark$expr == "BFGS"]
    time_BFGS_gr <- times[temp_microbenchmark$expr == "BFGS_gr"]
    time_L_BFGS_B <- times[temp_microbenchmark$expr == "L-BFGS-B"]
    time_L_BFGS_B_gr <- times[temp_microbenchmark$expr == "L-BFGS-B_gr"]
    time_CG <- times[temp_microbenchmark$expr == "CG"]
    time_CG_gr <- times[temp_microbenchmark$expr == "CG_gr"]
    time_Nelder_Mead <- times[temp_microbenchmark$expr == "Nelder-Mead"]
    time_nlm <- times[temp_microbenchmark$expr == "nlm"]
    time_nlm_gr <- times[temp_microbenchmark$expr == "nlm_gr"]
    time_nlm_he <- times[temp_microbenchmark$expr == "nlm_he"]
    time_nlm_grhe <- times[temp_microbenchmark$expr == "nlm_grhe"]
    time_nlminb <- times[temp_microbenchmark$expr == "nlminb"]
    time_nlminb_gr <- times[temp_microbenchmark$expr == "nlminb_gr"]
    time_nlminb_he <- times[temp_microbenchmark$expr == "nlminb_he"]
    time_nlminb_grhe <- times[temp_microbenchmark$expr == "nlminb_grhe"]
    time_hjn <- times[temp_microbenchmark$expr == "hjn"]
    # time_ucminf <- times[temp_microbenchmark$expr == "ucminf"]
    time_marqLevAlg <- times[temp_microbenchmark$expr == "marqLevAlg"]
    time_marqLevAlg_gr <- times[temp_microbenchmark$expr == "marqLevAlg_gr"]
    time_marqLevAlg_he <- times[temp_microbenchmark$expr == "marqLevAlg_he"]
    time_marqLevAlg_grhe <- times[temp_microbenchmark$expr == "marqLevAlg_grhe"]
    time_newuoa <- times[temp_microbenchmark$expr == "newuoa"]
    time_BBoptim <- times[temp_microbenchmark$expr == "BBoptim"]
    # time_BBoptim_gr <- times[temp_microbenchmark$expr == "BBoptim_gr"]
    # time_nloptr <- temp_microbenchmark_nloptr$time / 10^9
    
    # Get results -----------------
    result_BFGS_benchmark <- benchmark_model$all_natural_parameters_list$BFGS
    result_BFGS_gr_benchmark <- benchmark_model$all_natural_parameters_list$BFGS_gr
    result_L_BFGS_B_benchmark <- benchmark_model$all_natural_parameters_list$L_BFGS_B
    result_L_BFGS_B_gr_benchmark <- benchmark_model$all_natural_parameters_list$L_BFGS_B_gr
    result_CG_benchmark <- benchmark_model$all_natural_parameters_list$CG
    result_CG_gr_benchmark <- benchmark_model$all_natural_parameters_list$CG_gr
    result_Nelder_Mead_benchmark <- benchmark_model$all_natural_parameters_list$Nelder_Mead
    result_nlm_benchmark <- benchmark_model$all_natural_parameters_list$nlm
    result_nlm_gr_benchmark <- benchmark_model$all_natural_parameters_list$nlm_gr
    result_nlm_he_benchmark <- benchmark_model$all_natural_parameters_list$nlm_he
    result_nlm_grhe_benchmark <- benchmark_model$all_natural_parameters_list$nlm_grhe
    result_nlminb_benchmark <- benchmark_model$all_natural_parameters_list$nlminb
    result_nlminb_gr_benchmark <- benchmark_model$all_natural_parameters_list$nlminb_gr
    result_nlminb_he_benchmark <- benchmark_model$all_natural_parameters_list$nlminb_he
    result_nlminb_grhe_benchmark <- benchmark_model$all_natural_parameters_list$nlminb_grhe
    result_hjn_benchmark <- benchmark_model$all_natural_parameters_list$hjn
    result_marqLevAlg_benchmark <- benchmark_model$all_natural_parameters_list$marqLevAlg
    result_marqLevAlg_gr_benchmark <- benchmark_model$all_natural_parameters_list$marqLevAlg_gr
    result_marqLevAlg_he_benchmark <- benchmark_model$all_natural_parameters_list$marqLevAlg_he
    result_marqLevAlg_grhe_benchmark <- benchmark_model$all_natural_parameters_list$marqLevAlg_grhe
    result_newuoa_benchmark <- benchmark_model$all_natural_parameters_list$newuoa
    result_BBoptim_benchmark <- benchmark_model$all_natural_parameters_list$BBoptim
    
    # Store times and results ---------------------
    convergence_problems <- as_tibble(rbind(benchmark_model$problems)) %>%
      select(contains("failures")) %>%
      as.numeric()
    null_problems <- as_tibble(rbind(benchmark_model$problems)) %>%
      select(contains("null")) %>%
      as.numeric()
    
    params_mu <- setNames(lapply(1:m, function(x) bquote(c(result_BFGS_benchmark$mu[.(x)],
                                                           result_BFGS_gr_benchmark$mu[.(x)],
                                                           result_L_BFGS_B_benchmark$mu[.(x)],
                                                           result_L_BFGS_B_gr_benchmark$mu[.(x)],
                                                           result_CG_benchmark$mu[.(x)],
                                                           result_CG_gr_benchmark$mu[.(x)],
                                                           result_Nelder_Mead_benchmark$mu[.(x)],
                                                           result_nlm_benchmark$mu[.(x)],
                                                           result_nlm_gr_benchmark$mu[.(x)],
                                                           result_nlm_he_benchmark$mu[.(x)],
                                                           result_nlm_grhe_benchmark$mu[.(x)],
                                                           result_nlminb_benchmark$mu[.(x)],
                                                           result_nlminb_gr_benchmark$mu[.(x)],
                                                           result_nlminb_he_benchmark$mu[.(x)],
                                                           result_nlminb_grhe_benchmark$mu[.(x)],
                                                           result_hjn_benchmark$mu[.(x)],
                                                           result_marqLevAlg_benchmark$mu[.(x)],
                                                           result_marqLevAlg_gr_benchmark$mu[.(x)],
                                                           result_marqLevAlg_he_benchmark$mu[.(x)],
                                                           result_marqLevAlg_grhe_benchmark$mu[.(x)],
                                                           # result_ucminf_benchmark$mu[.(x)],
                                                           result_newuoa_benchmark$mu[.(x)],
                                                           result_BBoptim_benchmark$mu[.(x)]))),
                          mu_names)
    
    params_sigma <- setNames(lapply(1:m, function(x) bquote(c(result_BFGS_benchmark$sigma[.(x)],
                                                              result_BFGS_gr_benchmark$sigma[.(x)],
                                                              result_L_BFGS_B_benchmark$sigma[.(x)],
                                                              result_L_BFGS_B_gr_benchmark$sigma[.(x)],
                                                              result_CG_benchmark$sigma[.(x)],
                                                              result_CG_gr_benchmark$sigma[.(x)],
                                                              result_Nelder_Mead_benchmark$sigma[.(x)],
                                                              result_nlm_benchmark$sigma[.(x)],
                                                              result_nlm_gr_benchmark$sigma[.(x)],
                                                              result_nlm_he_benchmark$sigma[.(x)],
                                                              result_nlm_grhe_benchmark$sigma[.(x)],
                                                              result_nlminb_benchmark$sigma[.(x)],
                                                              result_nlminb_gr_benchmark$sigma[.(x)],
                                                              result_nlminb_he_benchmark$sigma[.(x)],
                                                              result_nlminb_grhe_benchmark$sigma[.(x)],
                                                              result_hjn_benchmark$sigma[.(x)],
                                                              result_marqLevAlg_benchmark$sigma[.(x)],
                                                              result_marqLevAlg_gr_benchmark$sigma[.(x)],
                                                              result_marqLevAlg_he_benchmark$sigma[.(x)],
                                                              result_marqLevAlg_grhe_benchmark$sigma[.(x)],
                                                              # result_ucminf_benchmark$sigma[.(x)],
                                                              result_newuoa_benchmark$sigma[.(x)],
                                                              result_BBoptim_benchmark$sigma[.(x)]))),
                             sigma_names)
    
    # params_lambda <- setNames(lapply(1:m, function(x) bquote(sapply(lapply(par_benchmark_list, "[[", "lambda"), "[", .(x)))),
    #                           accuracy_lambda_names)
    params_gamma <- setNames(lapply(1:(m^2), function(x) bquote(c(t(result_BFGS_benchmark$gamma)[.(x)],
                                                                  t(result_BFGS_gr_benchmark$gamma)[.(x)],
                                                                  t(result_L_BFGS_B_benchmark$gamma)[.(x)],
                                                                  t(result_L_BFGS_B_gr_benchmark$gamma)[.(x)],
                                                                  t(result_CG_benchmark$gamma)[.(x)],
                                                                  t(result_CG_gr_benchmark$gamma)[.(x)],
                                                                  t(result_Nelder_Mead_benchmark$gamma)[.(x)],
                                                                  t(result_nlm_benchmark$gamma)[.(x)],
                                                                  t(result_nlm_gr_benchmark$gamma)[.(x)],
                                                                  t(result_nlm_he_benchmark$gamma)[.(x)],
                                                                  t(result_nlm_grhe_benchmark$gamma)[.(x)],
                                                                  t(result_nlminb_benchmark$gamma)[.(x)],
                                                                  t(result_nlminb_gr_benchmark$gamma)[.(x)],
                                                                  t(result_nlminb_he_benchmark$gamma)[.(x)],
                                                                  t(result_nlminb_grhe_benchmark$gamma)[.(x)],
                                                                  t(result_hjn_benchmark$gamma)[.(x)],
                                                                  t(result_marqLevAlg_benchmark$gamma)[.(x)],
                                                                  t(result_marqLevAlg_gr_benchmark$gamma)[.(x)],
                                                                  t(result_marqLevAlg_he_benchmark$gamma)[.(x)],
                                                                  t(result_marqLevAlg_grhe_benchmark$gamma)[.(x)],
                                                                  # t(result_ucminf_benchmark$gamma)[.(x)],
                                                                  t(result_newuoa_benchmark$gamma)[.(x)],
                                                                  t(result_BBoptim_benchmark$gamma)[.(x)]))),
                             gamma_names)
    
    # params_gamma <- setNames(lapply(1:(m^2), function(x) bquote(sapply(lapply(lapply(par_benchmark_list, "[[", "gamma"), t), "[", .(x)))),
    #                          gamma_names)
    benchmark_optimizer_comparison_df_simu3 <- rbind(benchmark_optimizer_comparison_df_simu3,
                                                     cbind(data.frame(time = c(time_BFGS,
                                                                               time_BFGS_gr,
                                                                               time_L_BFGS_B,
                                                                               time_L_BFGS_B_gr,
                                                                               time_CG,
                                                                               time_CG_gr,
                                                                               time_Nelder_Mead,
                                                                               time_nlm,
                                                                               time_nlm_gr,
                                                                               time_nlm_he,
                                                                               time_nlm_grhe,
                                                                               time_nlminb,
                                                                               time_nlminb_gr,
                                                                               time_nlminb_he,
                                                                               time_nlminb_grhe,
                                                                               time_hjn,
                                                                               time_marqLevAlg,
                                                                               time_marqLevAlg_gr,
                                                                               time_marqLevAlg_he,
                                                                               time_marqLevAlg_grhe,
                                                                               # time_ucminf,
                                                                               time_newuoa,
                                                                               time_BBoptim),
                                                                      # time_nloptr),
                                                                      m = rep(m,
                                                                              # length(OPTIMIZERS_METHOD) + length(temp_microbenchmark_nloptr$expr)),
                                                                              length(OPTIMIZERS_METHOD)),
                                                                      optimizer = OPTIMIZERS_METHOD,
                                                                      # optimizer = c(OPTIMIZERS_METHOD, as.character(temp_microbenchmark_nloptr$expr)),
                                                                      iterations = c(benchmark_model$all_iterations$BFGS,
                                                                                     benchmark_model$all_iterations$BFGS_gr,
                                                                                     benchmark_model$all_iterations$L_BFGS_B,
                                                                                     benchmark_model$all_iterations$L_BFGS_B_gr,
                                                                                     benchmark_model$all_iterations$CG,
                                                                                     benchmark_model$all_iterations$CG_gr,
                                                                                     benchmark_model$all_iterations$Nelder_Mead,
                                                                                     benchmark_model$all_iterations$nlm,
                                                                                     benchmark_model$all_iterations$nlm_gr,
                                                                                     benchmark_model$all_iterations$nlm_he,
                                                                                     benchmark_model$all_iterations$nlm_grhe,
                                                                                     benchmark_model$all_iterations$nlminb,
                                                                                     benchmark_model$all_iterations$nlminb_gr,
                                                                                     benchmark_model$all_iterations$nlminb_he,
                                                                                     benchmark_model$all_iterations$nlminb_grhe,
                                                                                     benchmark_model$all_iterations$hjn,
                                                                                     benchmark_model$all_iterations$marqLevAlg,
                                                                                     benchmark_model$all_iterations$marqLevAlg_gr,
                                                                                     benchmark_model$all_iterations$marqLevAlg_he,
                                                                                     benchmark_model$all_iterations$marqLevAlg_grhe,
                                                                                     # benchmark_model$all_iterations$ucminf,
                                                                                     benchmark_model$all_iterations$newuoa,
                                                                                     benchmark_model$all_iterations$BBoptim),
                                                                      # result_nloptr_benchmark$iterations),
                                                                      nll = c(benchmark_model$all_nlls$BFGS,
                                                                              benchmark_model$all_nlls$BFGS_gr,
                                                                              benchmark_model$all_nlls$L_BFGS_B,
                                                                              benchmark_model$all_nlls$L_BFGS_B_gr,
                                                                              benchmark_model$all_nlls$CG,
                                                                              benchmark_model$all_nlls$CG_gr,
                                                                              benchmark_model$all_nlls$Nelder_Mead,
                                                                              benchmark_model$all_nlls$nlm,
                                                                              benchmark_model$all_nlls$nlm_gr,
                                                                              benchmark_model$all_nlls$nlm_he,
                                                                              benchmark_model$all_nlls$nlm_grhe,
                                                                              benchmark_model$all_nlls$nlminb,
                                                                              benchmark_model$all_nlls$nlminb_gr,
                                                                              benchmark_model$all_nlls$nlminb_he,
                                                                              benchmark_model$all_nlls$nlminb_grhe,
                                                                              benchmark_model$all_nlls$hjn,
                                                                              benchmark_model$all_nlls$marqLevAlg,
                                                                              benchmark_model$all_nlls$marqLevAlg_gr,
                                                                              benchmark_model$all_nlls$marqLevAlg_he,
                                                                              benchmark_model$all_nlls$marqLevAlg_grhe,
                                                                              # benchmark_model$all_nlls$ucminf,
                                                                              benchmark_model$all_nlls$newuoa,
                                                                              benchmark_model$all_nlls$BBoptim),
                                                                      # result_nloptr_benchmark$nll),
                                                                      convergence_problems = convergence_problems,
                                                                      null_problems = null_problems,
                                                                      seed = idx_counter_benchmark),
                                                           do.call("data.frame", append(params_mu, append(params_sigma, params_gamma)))))
    
    end_benchmark <- Sys.time()
    percent <- idx_counter_benchmark / BENCHMARK_SAMPLES
    if (percent %in% seq(0, 1, by = 0.1)) {
      duration <- end_benchmark - begin_benchmark
      total_time_to_completion <- duration / percent
      notif(title = "Simu3 Benchmark",
            paste0(100 * max(percent),
                   "%\ntotal time required: ", round(total_time_to_completion, 1), " ", units(total_time_to_completion),
                   "\nm: ", m,
                   "\nSize: ", DATA_SIZE_SIMU3))
    }
  }
  end_benchmark <- Sys.time()
  # percent <- BOOTSTRAP_SAMPLES / 1000
  duration <- end_benchmark - begin_benchmark
  # total_time_to_completion <- duration / percent
  
  notif(title = "Simu3 Benchmark",
        paste0(# 100 * percent,
          # "%\ntotal time required bootstrap: ", round(total_time_to_completion, 1), " ", units(total_time_to_completion),
          "duration: ", round(duration, 1), " ", units(duration),
          "\nm: ", m,
          "\nSize: ", DATA_SIZE_SIMU3))
}

# Bootstrap ---------------------------
begin_bootstrap <- Sys.time()
if (BOOTSTRAP_SAMPLES != 0) {
  # registerDoParallel(cores = CORES)
  registerDoFuture()
  plan(tweak(multisession, workers = CORES))
  bootstrap_simu3_results <- foreach(idx_counter_bootstrap = 1:BOOTSTRAP_SAMPLES,
                                     .packages = OPTIMIZER_PACKAGES,
                                     .combine = 'combine_foreach_rbind',
                                     .multicombine = TRUE) %dorng% {
                                       
                                       set.seed(idx_counter_bootstrap)
                                       # dyn.load(dynlib("code/poi_hmm_MSE_estimators"))
                                       dyn.load(dynlib("code/norm_hmm"))
                                       
                                       bootstrap_model <- norm.HMM.generate.estimable.sample(
                                         ns = DATA_SIZE_SIMU3,
                                         mod = list(m = m,
                                                    mu = simu3_true_mu,
                                                    sigma = simu3_true_sigma,
                                                    gamma = simu3_true_gamma,
                                                    delta = simu3_true_delta),
                                         # return_the_MSE_estimators = TRUE,
                                         debug_message = "bootstrap simu3 parallel"
                                       )
                                       
                                       # The values from gamma are taken columnwise
                                       natural_parameters <- bootstrap_model$natural_parameters
                                       natural_parameters <- unlist(natural_parameters[PARAMS_NAMES_NORM])
                                       natural_parameters <- c(natural_parameters, seed = idx_counter_bootstrap)
                                       all_natural_parameters <- bootstrap_model$all_natural_parameters
                                       all_natural_parameters <- cbind(all_natural_parameters, seed = idx_counter_bootstrap)
                                       all_nlls <- bootstrap_model$all_nlls
                                       all_nlls <- cbind(all_nlls, seed = idx_counter_bootstrap)
                                       problems <- bootstrap_model$problems
                                       problems <- c(problems, seed = idx_counter_bootstrap)
                                       # MSE_estimators <- bootstrap_model$all_MSE_estimators
                                       # MSE_estimators <- cbind(MSE_estimators, seed = idx_counter_bootstrap)
                                       
                                       # dyn.unload(dynlib("code/poi_hmm_MSE_estimators"))
                                       dyn.unload(dynlib("code/norm_hmm"))
                                       return(list(natural_parameters = natural_parameters,
                                                   all_natural_parameters = all_natural_parameters,
                                                   all_nlls = all_nlls,
                                                   problems = problems))
                                       # MSE_estimators = MSE_estimators))
                                     }
  # stopImplicitCluster()
  future:::ClusterRegistry("stop")
  # if (BOOTSTRAP_SAMPLES == 1) {
  #   names(bootstrap_simu3) <- params_names_latex
  # } else {
  #   colnames(bootstrap_simu3) <- params_names_latex
  # }
  # q <- apply(bootstrap_simu3,
  #            2,
  #            quantile.colwise)
  # conf_int_simu3$Bootstrap.L <- q[1, ]
  # conf_int_simu3$Bootstrap.U <- q[2, ]
}
end_bootstrap <- Sys.time()
percent <- BOOTSTRAP_SAMPLES / 1000
duration <- end_bootstrap - begin_bootstrap
total_time_to_completion <- duration / percent

notif(title = "Simu3 Bootstrap",
      paste0(100 * percent,
             "%\ntotal time required bootstrap: ", round(total_time_to_completion, 1), " ", units(total_time_to_completion),
             "\nm: ", m,
             "\nSize: ", DATA_SIZE_SIMU3))


# Robustness ---------------------------
set.seed(SEED)

# Replicate parameter_possibilities for all the simulations' true parameter settings
# New (dataset + true nlls) for each set of parameters
simu3_robustness <- lapply(X = true_params_simu3,
                           FUN = function(x, data_size, m) {
                             simu3_true_mu <- x$mu
                             simu3_true_sigma <- x$sigma
                             simu3_true_gamma <- x$gamma
                             simu3_true_delta <- x$delta
                             temp <- norm.HMM.generate.estimable.sample(ns = data_size,
                                                                        mod = list(m = m,
                                                                                   mu = simu3_true_mu,
                                                                                   sigma = simu3_true_sigma,
                                                                                   gamma = simu3_true_gamma,
                                                                                   delta = simu3_true_delta))
                             x$data <- temp$data
                             x$nlls <- temp$all_nlls
                             return(x)
                           },
                           data_size = DATA_SIZE_SIMU3,
                           m = M_LIST_SIMU3)

# For each new simulated data set, we need the parameter possibilities
unlink("mylogsimu3.txt")
registerDoFuture()
plan(multicore, workers = CORES)
parameter_possibilities_per_data_set <- foreach(dataset_idx = seq_along(true_params_simu3),
                                                .packages = packages,
                                                .combine = "combine_foreach_rbind",
                                                # .combine = "rbind",
                                                .inorder = TRUE) %dopar% {
                                                  
                                                  set.seed(SEED)
                                                  
                                                  simu_data <- simu3_robustness[[dataset_idx]]$data
                                                  
                                                  # Create parameters variations -------------------
                                                  accuracy_mus <- seq(min(simu_data), max(simu_data), by = 0.5)
                                                  max_temp <- max(simu_data)
                                                  mean_temp <- mean(simu_data)
                                                  min_temp <- min(simu_data)
                                                  # Bhatia-Davis inequality
                                                  upper_bound <- sqrt((max_temp - mean_temp) * (mean_temp - min_temp))
                                                  # von Szokefalvi Nagy inequality
                                                  lower_bound <- sqrt((max_temp - min_temp)^2 / (2 * length(simu_data)))
                                                  accuracy_sigmas <- seq(lower_bound, upper_bound, length.out = 10)
                                                  accuracy_gammas <- (1:9) / 10
                                                  
                                                  # params_lambda <- setNames(lapply(seq(-(m - 1), 0), function(x) bquote(head.secure(accuracy_lambdas, .(x)))), # lambda_1 stops m-1 steps before the end
                                                  #                           accuracy_lambda_names)
                                                  params_mu <- setNames(lapply(1:m, function(x) bquote(accuracy_mus)),
                                                                        mu_names)
                                                  params_sigma <- setNames(lapply(1:m, function(x) bquote(accuracy_sigmas)),
                                                                           sigma_names)
                                                  params_gamma <- setNames(lapply(1:(m^2), function(x) bquote(accuracy_gammas)),
                                                                           gamma_names)
                                                  
                                                  # Generate possible mus and sigmas -------------------
                                                  parameter_possibilities_mu <- do.call("CJ", params_mu)
                                                  parameter_possibilities_sigma <- do.call("CJ", params_sigma)
                                                  
                                                  # THEY ARE ALREADY IN INCREASING ORDER
                                                  # # Remove rows where lambdas are not in increasing order
                                                  # rows_to_keep <- t(apply(X = parameter_possibilities[, 1:m],
                                                  #                         MARGIN = 1,
                                                  #                         FUN = function(row) {
                                                  #                           m <- length(row)
                                                  #                           all(order(row) == 1:m & !any(duplicated(row))) # Keep only strictly increasing rows (delete rows which contain duplicates i.e. non-strictly increasing)
                                                  #                         }))
                                                  # rows_to_keep <- which(rows_to_keep)
                                                  # parameter_possibilities <- parameter_possibilities[rows_to_keep, ]
                                                  
                                                  parameter_possibilities_mu$idx <- 1:nrow(parameter_possibilities_mu)
                                                  
                                                  # Remove rows where the mus are not strictly increasing
                                                  # This could be applied instead on the sigmas, but we arbitrarily pick the mus instead
                                                  idx_rows_to_remove <- apply(X = parameter_possibilities_mu,
                                                                              MARGIN = 1,
                                                                              FUN = function(row, m2 = m) {
                                                                                
                                                                                row_idx <- row[["idx"]]
                                                                                
                                                                                row_values_lambda <- row[names(row) != "idx"]
                                                                                # Are any mus not sorted (not increasing), and are any mus not unique (not strict order)
                                                                                if (any(duplicated(row_values_lambda)) == TRUE || any(row_values_lambda != sort(row_values_lambda))) {
                                                                                  return(row_idx)
                                                                                }
                                                                                
                                                                              })
                                                  # Remove NULL elements and convert to vector at the same time (vectors cannot contain NULL)
                                                  idx_rows_to_remove <- unlist(idx_rows_to_remove)
                                                  
                                                  parameter_possibilities_mu <- parameter_possibilities_mu[-idx_rows_to_remove, ]
                                                  
                                                  # Generate possible gammas ---------------------
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
                                                  
                                                  # Prepare possibilities --------------------------
                                                  # data.frame objects are easier to manipulate than data.table objects
                                                  parameter_possibilities_mu <- as.data.frame(parameter_possibilities_mu)
                                                  parameter_possibilities_sigma <- as.data.frame(parameter_possibilities_sigma)
                                                  parameter_possibilities_gamma <- as.data.frame(parameter_possibilities_gamma)
                                                  # Remove idx column
                                                  parameter_possibilities_mu <- parameter_possibilities_mu[, names(parameter_possibilities_mu) != "idx"]
                                                  # parameter_possibilities_sigma did not get an idx column
                                                  parameter_possibilities_gamma <- parameter_possibilities_gamma[, names(parameter_possibilities_gamma) != "idx"]
                                                  
                                                  
                                                  # # Repeat the parameter_possibilities for each robustness dataset. dataset_idx links each set of parameter_possibilities to a dataset.
                                                  # full_parameter_possibilities <- rbindlist(lapply(1:length(simu3_robustness), function(x) cbind(parameter_possibilities, "dataset_idx" = x)))
                                                  # # data.table causes issues when selecting with variables
                                                  # full_parameter_possibilities <- as.data.frame(full_parameter_possibilities)
                                                  
                                                  # parameter_possibilities$idx <- 1:nrow(parameter_possibilities)
                                                  # parameter_possibilities$dataset_idx <- dataset_idx
                                                  
                                                  
                                                  
                                                  # K-Means (to remove) -------------------------
                                                  # # (we cannot apply k-means on the TPM because the mean would likely not respect the
                                                  # # constraint of rows summing to 1, so we apply k-means on the working - unconstrained - TPM instead)
                                                  # # TPM_cols <- grepl("gamma", names(parameter_possibilities_gamma))
                                                  # tgamma_possibilities <- t(apply(X = parameter_possibilities_gamma,
                                                  #                                 MARGIN = 1,
                                                  #                                 FUN = function(row) {
                                                  #                                   gamma <- matrix(as.numeric(row), nrow = m, ncol = m, byrow = TRUE)
                                                  #                                   result <- gamma.n2w(m, gamma)
                                                  #                                   return(result)
                                                  #                                 }))
                                                  # 
                                                  # parameter_possibilities_gamma_tgamma <- cbind(parameter_possibilities_gamma, tgamma_possibilities) %>%
                                                  #   # Default names are V1, V2, ... on Linux
                                                  #   rename_with(~ gsub("V", "tgamma", .x, fixed = TRUE), starts_with("V")) %>%
                                                  #   # Default names are 1, 2, ... on Windows
                                                  #   rename_with(~ gsub("^", "tgamma", .x), num_range(prefix = "", range = 1:(m^2 - m)))
                                                  # 
                                                  # # Cross join (i.e. cartesian join)
                                                  # parameter_possibilities <- merge(x = parameter_possibilities_mu, y = parameter_possibilities_sigma)
                                                  # parameter_possibilities <- merge(x = parameter_possibilities, y = parameter_possibilities_gamma_tgamma)
                                                  # 
                                                  # 
                                                  # # The default number of maximum iterations is too low. 1000 max iterations seems large enough.
                                                  # # The default algorithm (Hartigan-Wong) makes k-means throw warnings, so we use the second algorithm available in the list (Lloyd, see kmeans help)
                                                  # # 400 000 data in total should be manageable. With 16 data sets, each one needs to have a size of 400 000 / 16 = 25 000
                                                  # parameter_possibilities_means <- parameter_possibilities %>%
                                                  #   select(starts_with(c("mu", "sigma", "tgamma"))) %>%
                                                  #   kmeans(centers = 12500,
                                                  #          # algorithm = "Lloyd",
                                                  #          iter.max = 1000)
                                                  # # Mclust(G = 5)
                                                  # 
                                                  # # parameter_possibilities takes a lot of space (777 MB in one instance)
                                                  # rm(parameter_possibilities)
                                                  # 
                                                  # # Working tgamma -> natural gamma
                                                  # result <- foreach(i = 1:nrow(parameter_possibilities_means$centers),
                                                  #                   .combine = "rbind",
                                                  #                   .inorder = TRUE) %do% {
                                                  #                     
                                                  #                     gamma <- gamma.w2n(m, parameter_possibilities_means$centers[i, tgamma_indices])
                                                  #                     gamma_vector <- as.vector(t(gamma)) # row-wise instead of default column-wise
                                                  #                     mu <- parameter_possibilities_means$centers[i, mu_indices]
                                                  #                     sigma <- parameter_possibilities_means$centers[i, sigma_indices]
                                                  #                     result <- c(mu, sigma, gamma_vector, dataset_idx = dataset_idx)
                                                  #                     return(result)
                                                  #                   }
                                                  # 
                                                  # return(append(list(result = result), parameter_possibilities_means))
                                                  
                                                  # Sampling -------------------
                                                  # Cross join (i.e. cartesian join)
                                                  parameter_possibilities <- merge(x = parameter_possibilities_mu, y = parameter_possibilities_sigma)
                                                  parameter_possibilities <- merge(x = parameter_possibilities, y = parameter_possibilities_gamma)
                                                  parameter_possibilities_initial_amounts <- nrow(parameter_possibilities)
                                                  # Sample 12 500 initial values (we want 200 000 values in total, about the same number as in other simulations)
                                                  # 200 000 / 16 = 12 500
                                                  sample_indices <- sample(x = 1:nrow(parameter_possibilities), size = 12500)
                                                  parameter_possibilities <- parameter_possibilities[sample_indices, ]
                                                  parameter_possibilities$dataset_idx <- dataset_idx
                                                  
                                                  return(list(parameter_possibilities = parameter_possibilities, parameter_possibilities_initial_amounts = parameter_possibilities_initial_amounts))
                                                }
future:::ClusterRegistry("stop")

len_temp <- length(mu_names) + length(sigma_names) + length(gamma_names)
# colnames(full_parameter_possibilities$result)[1:len_temp] <- c(mu_names, sigma_names, gamma_names)
# rownames(full_parameter_possibilities$result) <- NULL
# full_parameter_possibilities$result <- as.data.frame(full_parameter_possibilities$result)
# full_parameter_possibilities$result$idx <- 1:nrow(full_parameter_possibilities$result)
full_parameter_possibilities <- parameter_possibilities_per_data_set$parameter_possibilities
simu3_parameter_possibilities_initial_amounts <- parameter_possibilities_per_data_set$parameter_possibilities_initial_amount
colnames(full_parameter_possibilities)[1:len_temp] <- c(mu_names, sigma_names, gamma_names)
rownames(full_parameter_possibilities) <- NULL
full_parameter_possibilities <- as.data.frame(full_parameter_possibilities)
full_parameter_possibilities$idx <- 1:nrow(full_parameter_possibilities)

# simu3_kmeans_objects <- full_parameter_possibilities
# simu3_kmeans_objects$result <- NULL


begin_accuracy <- Sys.time()

# https://stackoverflow.com/a/33632696/14323031
# cl <- makePSOCKcluster(CORES, outfile = "")
# registerDoParallel(cl)
# registerDoParallel(cores = CORES)
# pb <- txtProgressBar(1, max(full_parameter_possibilities$idx), style = 3)
full_parameter_possibilities_max_idx <- max(full_parameter_possibilities$idx)

# unlink("mylogsimu3.txt")
# file.create("mylogsimu3.txt")
# lockfile_path <- paste0(getwd(), "/lock")

registerDoFuture()
if (Sys.info()[['sysname']] == "Windows") {
  plan(multisession, workers = CORES)
} else {
  plan(multicore, workers = CORES)
}
# Some tasks (processes) will finish sooner than others. If processes are idle, then the whole task takes longer.
# So we split the main task into smaller tasks.
# 4 * 32 cores = 128 tasks. This should be enough to avoid too much idleness, but not too large so that
# too much time is wasted on pre-processing such as forking and loading packages.
# There are approximately 100 000 or more rows to process. 100k/128 =~ 800 rows to work on per process.
# results <- foreach(parameter_possibilities_rows = isplitRows(full_parameter_possibilities[1:10, ], chunks = 4 * CORES),
results <- foreach(parameter_possibilities_row = iter(full_parameter_possibilities, by = "row"),
                   .packages = c(OPTIMIZER_PACKAGES, "flock"),
                   .options.future = list(chunk.size = 15),
                   .combine = rbind,
                   # .combine = 'combine_foreach_rbind',
                   # .multicombine = TRUE,
                   .inorder = TRUE) %dorng% {
                     
                     # lock_log <- lock(lockfile_path)
                     # cat(paste0(parameter_possibilities_row$idx, "\n"), sep = "", file="mylogsimu3.txt", append=TRUE)
                     # unlock(lock_log)
                     
                     # setTxtProgressBar(pb, max(parameter_possibilities_rows$idx))
                     set.seed(SEED)
                     dyn.load(dynlib("code/norm_hmm"))
                     
                     # Variables ---------
                     # accuracy_samples_amount <- nrow(parameter_possibilities_rows)
                     
                     accuracy_rates_names <- c(paste0(OPTIMIZERS_METHOD_UNDERSCORE, "_convergence_failure"), paste0(OPTIMIZERS_METHOD_UNDERSCORE, "_nll_found"))
                     # accuracy_rates_names <- setNames(lapply(accuracy_rates_names, function(x) bquote(rep(NA, accuracy_samples_amount))),
                     #                                  accuracy_rates_names)
                     accuracy_rates_names <- setNames(lapply(accuracy_rates_names, function(x) NA),
                                                      accuracy_rates_names)
                     accuracy_rates_simu3 <- do.call(data.frame, accuracy_rates_names)
                     accuracy_rates_simu3 <- accuracy_rates_simu3[, order(names(accuracy_rates_simu3))]
                     # accuracy_rates_simu3 <- cbind(m = rep(m, accuracy_samples_amount),
                     #                              # dataset_idx = parameter_possibilities_rows$dataset_idx,
                     #                              accuracy_rates_simu3)
                     accuracy_rates_simu3 <- cbind(m = m,
                                                   dataset_idx = parameter_possibilities_row$dataset_idx,
                                                   accuracy_rates_simu3,
                                                   idx = parameter_possibilities_row$idx)
                     
                     # accuracy_mles_simu3 <- data.frame()
                     
                     accuracy_dataset_idx <- parameter_possibilities_row$dataset_idx
                     accuracy_variables <- simu3_robustness[[accuracy_dataset_idx]]
                     
                     accuracy_dataset <- accuracy_variables$data
                     reference_nlls_accuracy <- accuracy_variables$nlls
                     # There are some rare instances of abnormally high nlls. The median is a good way to get a "true" nll
                     reference_nll_accuracy <- median(as.numeric(reference_nlls_accuracy))
                     reference_mu_accuracy <- accuracy_variables$mu
                     reference_sigma_accuracy <- accuracy_variables$sigma
                     reference_gamma_accuracy <- accuracy_variables$gamma
                     reference_delta_accuracy <- accuracy_variables$delta
                     
                     # Starting parameters with some variation -----------
                     accuracy_mu <- as.numeric(parameter_possibilities_row[, mu_indices])
                     accuracy_sigma <- as.numeric(parameter_possibilities_row[, sigma_indices])
                     accuracy_gamma_temp <- as.numeric(parameter_possibilities_row[, gamma_indices])
                     accuracy_gamma <- matrix(accuracy_gamma_temp, nrow = m, ncol = m, byrow = TRUE)
                     
                     # Parameters & covariates for DM and TMB ----------
                     accuracy_working_params <- norm.HMM.pn2pw(m, accuracy_mu, accuracy_sigma, accuracy_gamma)
                     TMB_accuracy_data <- list(x = accuracy_dataset,
                                               m = m)
                     
                     # Optimization --------
                     # We can ignore all text
                     # sink(nullfile())
                     result_BFGS_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                               working_parameters = accuracy_working_params,
                                                               optimizer = "BFGS",
                                                               debug_message = "Simu3 > accuracy")
                     result_BFGS_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                  working_parameters = accuracy_working_params,
                                                                  gradient = TRUE,
                                                                  optimizer = "BFGS",
                                                                  debug_message = "Simu3 > accuracy")
                     result_L_BFGS_B_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                   working_parameters = accuracy_working_params,
                                                                   optimizer = "L-BFGS-B",
                                                                   debug_message = "Simu3 > accuracy")
                     result_L_BFGS_B_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                      working_parameters = accuracy_working_params,
                                                                      gradient = TRUE,
                                                                      optimizer = "L-BFGS-B",
                                                                      debug_message = "Simu3 > accuracy")
                     result_CG_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                             working_parameters = accuracy_working_params,
                                                             optimizer = "CG",
                                                             debug_message = "Simu3 > accuracy")
                     result_CG_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                working_parameters = accuracy_working_params,
                                                                gradient = TRUE,
                                                                optimizer = "CG",
                                                                debug_message = "Simu3 > accuracy")
                     result_Nelder_Mead_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                      working_parameters = accuracy_working_params,
                                                                      optimizer = "Nelder-Mead",
                                                                      debug_message = "Simu3 > accuracy")
                     result_nlm_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                              working_parameters = accuracy_working_params,
                                                              optimizer = "nlm",
                                                              debug_message = "Simu3 > accuracy")
                     result_nlm_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                 working_parameters = accuracy_working_params,
                                                                 gradient = TRUE,
                                                                 optimizer = "nlm",
                                                                 debug_message = "Simu3 > accuracy")
                     result_nlm_he_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                 working_parameters = accuracy_working_params,
                                                                 hessian = TRUE,
                                                                 optimizer = "nlm",
                                                                 debug_message = "Simu3 > accuracy")
                     result_nlm_grhe_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                   working_parameters = accuracy_working_params,
                                                                   gradient = TRUE,
                                                                   hessian = TRUE,
                                                                   optimizer = "nlm",
                                                                   debug_message = "Simu3 > accuracy")
                     result_nlminb_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                 working_parameters = accuracy_working_params,
                                                                 optimizer = "nlminb",
                                                                 debug_message = "Simu3 > accuracy")
                     result_nlminb_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                    working_parameters = accuracy_working_params,
                                                                    gradient = TRUE,
                                                                    optimizer = "nlminb",
                                                                    debug_message = "Simu3 > accuracy")
                     result_nlminb_he_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                    working_parameters = accuracy_working_params,
                                                                    hessian = TRUE,
                                                                    optimizer = "nlminb",
                                                                    debug_message = "Simu3 > accuracy")
                     result_nlminb_grhe_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                      working_parameters = accuracy_working_params,
                                                                      gradient = TRUE,
                                                                      hessian = TRUE,
                                                                      optimizer = "nlminb",
                                                                      debug_message = "Simu3 > accuracy")
                     result_hjn_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                              working_parameters = accuracy_working_params,
                                                              optimizer = "hjn",
                                                              debug_message = "Simu3 > accuracy")
                     result_marqLevAlg_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                     working_parameters = accuracy_working_params,
                                                                     optimizer = "marqLevAlg",
                                                                     debug_message = "Simu3 > accuracy")
                     result_marqLevAlg_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                        working_parameters = accuracy_working_params,
                                                                        gradient = TRUE,
                                                                        optimizer = "marqLevAlg",
                                                                        debug_message = "Simu3 > accuracy")
                     result_marqLevAlg_he_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                        working_parameters = accuracy_working_params,
                                                                        hessian = TRUE,
                                                                        optimizer = "marqLevAlg",
                                                                        debug_message = "Simu3 > accuracy")
                     result_marqLevAlg_grhe_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                          working_parameters = accuracy_working_params,
                                                                          gradient = TRUE,
                                                                          hessian = TRUE,
                                                                          optimizer = "marqLevAlg",
                                                                          debug_message = "Simu3 > accuracy")
                     # result_ucminf_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                     #                               working_parameters = accuracy_working_params,
                     #                               gradient = TRUE,
                     #                               hessian = TRUE,
                     #                               optimizer = "ucminf",
                     #                               debug_message = "Simu3 > accuracy")
                     result_newuoa_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                 working_parameters = accuracy_working_params,
                                                                 optimizer = "newuoa",
                                                                 debug_message = "Simu3 > accuracy")
                     result_BBoptim_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                                                                  working_parameters = accuracy_working_params,
                                                                  optimizer = "BBoptim",
                                                                  debug_message = "Simu3 > accuracy")
                     # result_BBoptim_gr_accuracy <- norm.TMB.estimate(TMB_data = TMB_accuracy_data,
                     #                                   working_parameters = accuracy_working_params,
                     #                                   gradient = TRUE,
                     #                                   optimizer = "BBoptim",
                     #                                   debug_message = "Simu3 > accuracy")
                     
                     # Check convergence --------
                     accuracy_rates_simu3$BFGS_convergence_failure <- is.null(result_BFGS_accuracy) || !result_BFGS_accuracy$convergence
                     accuracy_rates_simu3$BFGS_gr_convergence_failure <- is.null(result_BFGS_gr_accuracy) || !result_BFGS_gr_accuracy$convergence
                     accuracy_rates_simu3$L_BFGS_B_convergence_failure <- is.null(result_L_BFGS_B_accuracy) || !result_L_BFGS_B_accuracy$convergence
                     accuracy_rates_simu3$L_BFGS_B_gr_convergence_failure <- is.null(result_L_BFGS_B_gr_accuracy) || !result_L_BFGS_B_gr_accuracy$convergence
                     accuracy_rates_simu3$CG_convergence_failure <- is.null(result_CG_accuracy) || !result_CG_accuracy$convergence
                     accuracy_rates_simu3$CG_gr_convergence_failure <- is.null(result_CG_gr_accuracy) || !result_CG_gr_accuracy$convergence
                     accuracy_rates_simu3$Nelder_Mead_convergence_failure <- is.null(result_Nelder_Mead_accuracy) || !result_Nelder_Mead_accuracy$convergence
                     accuracy_rates_simu3$nlm_convergence_failure <- is.null(result_nlm_accuracy) || !result_nlm_accuracy$convergence
                     accuracy_rates_simu3$nlm_gr_convergence_failure <- is.null(result_nlm_gr_accuracy) || !result_nlm_gr_accuracy$convergence
                     accuracy_rates_simu3$nlm_he_convergence_failure <- is.null(result_nlm_he_accuracy) || !result_nlm_he_accuracy$convergence
                     accuracy_rates_simu3$nlm_grhe_convergence_failure <- is.null(result_nlm_grhe_accuracy) || !result_nlm_grhe_accuracy$convergence
                     accuracy_rates_simu3$nlminb_convergence_failure <- is.null(result_nlminb_accuracy) || !result_nlminb_accuracy$convergence
                     accuracy_rates_simu3$nlminb_gr_convergence_failure <- is.null(result_nlminb_gr_accuracy) || !result_nlminb_gr_accuracy$convergence
                     accuracy_rates_simu3$nlminb_he_convergence_failure <- is.null(result_nlminb_he_accuracy) || !result_nlminb_he_accuracy$convergence
                     accuracy_rates_simu3$nlminb_grhe_convergence_failure <- is.null(result_nlminb_grhe_accuracy) || !result_nlminb_grhe_accuracy$convergence
                     accuracy_rates_simu3$hjn_convergence_failure <- is.null(result_hjn_accuracy) || !result_hjn_accuracy$convergence
                     accuracy_rates_simu3$marqLevAlg_convergence_failure <- is.null(result_marqLevAlg_accuracy) || !result_marqLevAlg_accuracy$convergence
                     accuracy_rates_simu3$marqLevAlg_gr_convergence_failure <- is.null(result_marqLevAlg_gr_accuracy) || !result_marqLevAlg_gr_accuracy$convergence
                     accuracy_rates_simu3$marqLevAlg_he_convergence_failure <- is.null(result_marqLevAlg_he_accuracy) || !result_marqLevAlg_he_accuracy$convergence
                     accuracy_rates_simu3$marqLevAlg_grhe_convergence_failure <- is.null(result_marqLevAlg_grhe_accuracy) || !result_marqLevAlg_grhe_accuracy$convergence
                     accuracy_rates_simu3$newuoa_convergence_failure <- is.null(result_newuoa_accuracy) || !result_newuoa_accuracy$convergence
                     accuracy_rates_simu3$BBoptim_convergence_failure <- is.null(result_BBoptim_accuracy) || !result_BBoptim_accuracy$convergence
                     
                     # Check nll ------------
                     accuracy_rates_simu3$BFGS_nll_found <- if (accuracy_rates_simu3$BFGS_convergence_failure) NA else abs(result_BFGS_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$BFGS_gr_nll_found <- if (accuracy_rates_simu3$BFGS_gr_convergence_failure) NA else abs(result_BFGS_gr_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$L_BFGS_B_nll_found <- if (accuracy_rates_simu3$L_BFGS_B_convergence_failure) NA else abs(result_L_BFGS_B_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$L_BFGS_B_gr_nll_found <- if (accuracy_rates_simu3$L_BFGS_B_gr_convergence_failure) NA else abs(result_L_BFGS_B_gr_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$CG_nll_found <- if (accuracy_rates_simu3$CG_convergence_failure) NA else abs(result_CG_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$CG_gr_nll_found <- if (accuracy_rates_simu3$CG_gr_convergence_failure) NA else abs(result_CG_gr_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$Nelder_Mead_nll_found <- if (accuracy_rates_simu3$Nelder_Mead_convergence_failure) NA else abs(result_Nelder_Mead_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlm_nll_found <- if (accuracy_rates_simu3$nlm_convergence_failure) NA else abs(result_nlm_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlm_gr_nll_found <- if (accuracy_rates_simu3$nlm_gr_convergence_failure) NA else abs(result_nlm_gr_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlm_he_nll_found <- if (accuracy_rates_simu3$nlm_he_convergence_failure) NA else abs(result_nlm_he_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlm_grhe_nll_found <- if (accuracy_rates_simu3$nlm_grhe_convergence_failure) NA else abs(result_nlm_grhe_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlminb_nll_found <- if (accuracy_rates_simu3$nlminb_convergence_failure) NA else abs(result_nlminb_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlminb_gr_nll_found <- if (accuracy_rates_simu3$nlminb_gr_convergence_failure) NA else abs(result_nlminb_gr_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlminb_he_nll_found <- if (accuracy_rates_simu3$nlminb_he_convergence_failure) NA else abs(result_nlminb_he_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$nlminb_grhe_nll_found <- if (accuracy_rates_simu3$nlminb_grhe_convergence_failure) NA else abs(result_nlminb_grhe_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$hjn_nll_found <- if (accuracy_rates_simu3$hjn_convergence_failure) NA else abs(result_hjn_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$marqLevAlg_nll_found <- if (accuracy_rates_simu3$marqLevAlg_convergence_failure) NA else abs(result_marqLevAlg_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$marqLevAlg_gr_nll_found <- if (accuracy_rates_simu3$marqLevAlg_gr_convergence_failure) NA else abs(result_marqLevAlg_gr_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$marqLevAlg_he_nll_found <- if (accuracy_rates_simu3$marqLevAlg_he_convergence_failure) NA else abs(result_marqLevAlg_he_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$marqLevAlg_grhe_nll_found <- if (accuracy_rates_simu3$marqLevAlg_grhe_convergence_failure) NA else abs(result_marqLevAlg_grhe_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$newuoa_nll_found <- if (accuracy_rates_simu3$newuoa_convergence_failure) NA else abs(result_newuoa_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     accuracy_rates_simu3$BBoptim_nll_found <- if (accuracy_rates_simu3$BBoptim_convergence_failure) NA else abs(result_BBoptim_accuracy$nll / reference_nll_accuracy - 1) <= NLL_THRESHOLD
                     
                     
                     # # Record estimates & results ---------
                     # # par_BFGS_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "BFGS_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_BFGS_accuracy$par) else NA
                     # # par_L_BFGS_B_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "L_BFGS_B_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_L_BFGS_B_accuracy$par) else NA
                     # # par_CG_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "CG_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_CG_accuracy$par) else NA
                     # # par_Nelder_Mead_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "Nelder_Mead_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_Nelder_Mead_accuracy$par) else NA
                     # # par_nlm_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlm_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_nlm_accuracy$estimate) else NA
                     # # par_nlminb_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlminb_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_nlminb_accuracy$par) else NA
                     # # par_hjn_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "hjn_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_hjn_accuracy$par) else NA
                     # # par_marqLevAlg_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "marqLevAlg_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_marqLevAlg_accuracy$b) else NA
                     # # # par_ucminf_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "ucminf_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_ucminf_accuracy$par) else NA
                     # # par_newuoa_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "newuoa_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_newuoa_accuracy$par) else NA
                     # # par_BBoptim_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "BBoptim_convergence_failure"] == FALSE) pois.HMM.pw2pn(m, result_BBoptim_accuracy$par) else NA
                     # 
                     # # par_BFGS_accuracy <- if (!is.na(par_BFGS_accuracy)) pois.HMM.label.order(m = m, lambda = par_BFGS_accuracy$lambda, gamma = par_BFGS_accuracy$gamma, delta = par_BFGS_accuracy$delta) else NA
                     # # par_L_BFGS_B_accuracy <- if (!is.na(par_L_BFGS_B_accuracy)) pois.HMM.label.order(m = m, lambda = par_L_BFGS_B_accuracy$lambda, gamma = par_L_BFGS_B_accuracy$gamma, delta = par_L_BFGS_B_accuracy$delta) else NA
                     # # par_CG_accuracy <- if (!is.na(par_CG_accuracy)) pois.HMM.label.order(m = m, lambda = par_CG_accuracy$lambda, gamma = par_CG_accuracy$gamma, delta = par_CG_accuracy$delta) else NA
                     # # par_Nelder_Mead_accuracy <- if (!is.na(par_Nelder_Mead_accuracy)) pois.HMM.label.order(m = m, lambda = par_Nelder_Mead_accuracy$lambda, gamma = par_Nelder_Mead_accuracy$gamma, delta = par_Nelder_Mead_accuracy$delta) else NA
                     # # par_nlm_accuracy <- if (!is.na(par_nlm_accuracy)) pois.HMM.label.order(m = m, lambda = par_nlm_accuracy$lambda, gamma = par_nlm_accuracy$gamma, delta = par_nlm_accuracy$delta) else NA
                     # # par_nlminb_accuracy <- if (!is.na(par_nlminb_accuracy)) pois.HMM.label.order(m = m, lambda = par_nlminb_accuracy$lambda, gamma = par_nlminb_accuracy$gamma, delta = par_nlminb_accuracy$delta) else NA
                     # # par_hjn_accuracy <- if (!is.na(par_hjn_accuracy)) pois.HMM.label.order(m = m, lambda = par_hjn_accuracy$lambda, gamma = par_hjn_accuracy$gamma, delta = par_hjn_accuracy$delta) else NA
                     # # par_marqLevAlg_accuracy <- if (!is.na(par_marqLevAlg_accuracy)) pois.HMM.label.order(m = m, lambda = par_marqLevAlg_accuracy$lambda, gamma = par_marqLevAlg_accuracy$gamma, delta = par_marqLevAlg_accuracy$delta) else NA
                     # # # par_ucminf_accuracy <- if (!is.na(par_ucminf_accuracy)) pois.HMM.label.order(m = m, lambda = par_ucminf_accuracy$lambda, gamma = par_ucminf_accuracy$gamma, delta = par_ucminf_accuracy$delta) else NA
                     # # par_newuoa_accuracy <- if (!is.na(par_newuoa_accuracy)) pois.HMM.label.order(m = m, lambda = par_newuoa_accuracy$lambda, gamma = par_newuoa_accuracy$gamma, delta = par_newuoa_accuracy$delta) else NA
                     # # par_BBoptim_accuracy <- if (!is.na(par_BBoptim_accuracy)) pois.HMM.label.order(m = m, lambda = par_BBoptim_accuracy$lambda, gamma = par_BBoptim_accuracy$gamma, delta = par_BBoptim_accuracy$delta) else NA
                     # 
                     # # Pass the estimates if there was no convergence failure
                     # estimates_BFGS_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "BFGS_convergence_failure"]) NA else unlist(result_BFGS_accuracy[PARAMS_NAMES_NORM])
                     # estimates_BFGS_gr_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "BFGS_gr_convergence_failure"]) NA else unlist(result_BFGS_gr_accuracy[PARAMS_NAMES_NORM])
                     # estimates_L_BFGS_B_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "L_BFGS_B_convergence_failure"]) NA else unlist(result_L_BFGS_B_accuracy[PARAMS_NAMES_NORM])
                     # estimates_L_BFGS_B_gr_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "L_BFGS_B_gr_convergence_failure"]) NA else unlist(result_L_BFGS_B_gr_accuracy[PARAMS_NAMES_NORM])
                     # estimates_CG_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "CG_convergence_failure"]) NA else unlist(result_CG_accuracy[PARAMS_NAMES_NORM])
                     # estimates_CG_gr_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "CG_gr_convergence_failure"]) NA else unlist(result_CG_gr_accuracy[PARAMS_NAMES_NORM])
                     # estimates_Nelder_Mead_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "Nelder_Mead_convergence_failure"]) NA else unlist(result_Nelder_Mead_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlm_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlm_convergence_failure"]) NA else unlist(result_nlm_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlm_gr_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlm_gr_convergence_failure"]) NA else unlist(result_nlm_gr_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlm_he_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlm_he_convergence_failure"]) NA else unlist(result_nlm_he_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlm_grhe_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlm_grhe_convergence_failure"]) NA else unlist(result_nlm_grhe_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlminb_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlminb_convergence_failure"]) NA else unlist(result_nlminb_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlminb_gr_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlminb_gr_convergence_failure"]) NA else unlist(result_nlminb_gr_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlminb_he_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlminb_he_convergence_failure"]) NA else unlist(result_nlminb_he_accuracy[PARAMS_NAMES_NORM])
                     # estimates_nlminb_grhe_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "nlminb_grhe_convergence_failure"]) NA else unlist(result_nlminb_grhe_accuracy[PARAMS_NAMES_NORM])
                     # estimates_hjn_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "hjn_convergence_failure"]) NA else unlist(result_hjn_accuracy[PARAMS_NAMES_NORM])
                     # estimates_marqLevAlg_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "marqLevAlg_convergence_failure"]) NA else unlist(result_marqLevAlg_accuracy[PARAMS_NAMES_NORM])
                     # estimates_marqLevAlg_gr_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "marqLevAlg_gr_convergence_failure"]) NA else unlist(result_marqLevAlg_gr_accuracy[PARAMS_NAMES_NORM])
                     # estimates_marqLevAlg_he_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "marqLevAlg_he_convergence_failure"]) NA else unlist(result_marqLevAlg_he_accuracy[PARAMS_NAMES_NORM])
                     # estimates_marqLevAlg_grhe_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "marqLevAlg_grhe_convergence_failure"]) NA else unlist(result_marqLevAlg_grhe_accuracy[PARAMS_NAMES_NORM])
                     # estimates_newuoa_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "newuoa_convergence_failure"]) NA else unlist(result_newuoa_accuracy[PARAMS_NAMES_NORM])
                     # estimates_BBoptim_accuracy <- if (accuracy_rates_simu3[idx_counter_accuracy, "BBoptim_convergence_failure"]) NA else unlist(result_BBoptim_accuracy[PARAMS_NAMES_NORM])
                     # 
                     # rownames_mles_accuracy <- c(paste0("lambda", 1:m), paste0("gamma", 1:(m ^ 2)), paste0("delta", 1:m))
                     # accuracy_mles_simu3 <- rbind(accuracy_mles_simu3, data.frame("m" = m,
                     #                                                              "par_name" = rownames_mles_accuracy,
                     #                                                              "par_reference" = c(reference_lambda_accuracy, reference_gamma_accuracy, reference_delta_accuracy),
                     #                                                              "par_initial" = c(accuracy_lambda, accuracy_gamma, stat.dist(accuracy_gamma)),
                     #                                                              "BFGS" = estimates_BFGS_accuracy,
                     #                                                              "BFGS_gr" = estimates_BFGS_gr_accuracy,
                     #                                                              "L_BFGS_B" = estimates_L_BFGS_B_accuracy,
                     #                                                              "L_BFGS_B_gr" = estimates_L_BFGS_B_gr_accuracy,
                     #                                                              "CG" = estimates_CG_accuracy,
                     #                                                              "CG_gr" = estimates_CG_gr_accuracy,
                     #                                                              "Nelder_Mead" = estimates_Nelder_Mead_accuracy,
                     #                                                              "nlm" = estimates_nlm_accuracy,
                     #                                                              "nlm_gr" = estimates_nlm_gr_accuracy,
                     #                                                              "nlm_he" = estimates_nlm_he_accuracy,
                     #                                                              "nlm_grhe" = estimates_nlm_grhe_accuracy,
                     #                                                              "nlminb" = estimates_nlminb_accuracy,
                     #                                                              "nlminb_gr" = estimates_nlminb_gr_accuracy,
                     #                                                              "nlminb_he" = estimates_nlminb_he_accuracy,
                     #                                                              "nlminb_grhe" = estimates_nlminb_grhe_accuracy,
                     #                                                              "hjn" = estimates_hjn_accuracy,
                     #                                                              "marqLevAlg" = estimates_marqLevAlg_accuracy,
                     #                                                              "marqLevAlg_gr" = estimates_marqLevAlg_gr_accuracy,
                     #                                                              "marqLevAlg_he" = estimates_marqLevAlg_he_accuracy,
                     #                                                              "marqLevAlg_grhe" = estimates_marqLevAlg_grhe_accuracy,
                     #                                                              # "ucminf" = estimates_ucminf_accuracy,
                     #                                                              "newuoa" = estimates_newuoa_accuracy,
                     #                                                              "BBoptim" = estimates_BBoptim_accuracy,
                     #                                                              "parameter_possibility_idx" = parameter_possibilities_rows$idx[idx_counter_accuracy],
                     #                                                              "dataset_idx" = dataset_accuracy_idx),
                     #                              make.row.names = FALSE)
                     
                     # Return --------------
                     dyn.unload(dynlib("code/norm_hmm"))
                     
                     end_accuracy <- Sys.time()
                     
                     # lock_log <- lock(lockfile_path)
                     # # if (Sys.info()[['sysname']] == "Windows") {
                     # #   shell(paste0("sed -ri 's/^(", parameter_possibilities_row$idx ,")$/\\1 END/' mylogsimu3.txt"))
                     # # } else {
                     # #   system(paste0("sed -ri 's/(", parameter_possibilities_row$idx ,")$/\\1 END/' mylogsimu3.txt"))
                     # # }
                     # cat(paste0(parameter_possibilities_row$idx, "\n"), sep = "", file="mylogsimu3.txt", append=TRUE)
                     # 
                     # if (Sys.info()[['sysname']] == "Windows") {
                     #   # current_amount_done <- as.numeric(shell("grep -c 'END' mylogsimu3.txt", intern = TRUE, mustWork = NA))
                     #   current_amount_done <- as.numeric(shell("wc -l mylogsimu3.txt | cut -f 1 -d ' '", intern = TRUE, mustWork = NA))
                     # } else {
                     #   # current_amount_done <- as.numeric(system("wc -l mylogsimu3.txt | cut -f 1 -d ' '", intern = TRUE))
                     #   current_amount_done <- as.numeric(system("wc -l mylogsimu3.txt | cut -f 1 -d ' '", intern = TRUE))
                     # }
                     # unlock(lock_log)
                     
                     current_amount_done <- parameter_possibilities_row$idx
                     percent <- current_amount_done / full_parameter_possibilities_max_idx
                     notification_targets_indices <- (1:10) * (nrow(full_parameter_possibilities) %/% 10)
                     if (current_amount_done %in% notification_targets_indices) {
                       duration <- end_accuracy - begin_accuracy
                       total_time_to_completion <- duration / max(percent)
                       time_left_to_completion <- total_time_to_completion - duration
                       notif(title = "Simu3 Robustness",
                             paste0(100 * round(percent, 2),
                                    "%\ntotal time robustness: ", round(total_time_to_completion, 1), " ", units(total_time_to_completion),
                                    "\ntime left robustness: ", round(time_left_to_completion, 1), " ", units(time_left_to_completion),
                                    "\nm: ", m,
                                    "\nSize: ", DATA_SIZE_SIMU3,
                                    "\nTotal amount:", nrow(full_parameter_possibilities)))
                     }
                     
                     # return(list(accuracy_rates_simu3 = accuracy_rates_simu3, accuracy_mles_simu3 = accuracy_mles_simu3))
                     # return(list(accuracy_rates_simu3 = accuracy_rates_simu3))
                     return(accuracy_rates_simu3)
                   }
# stopCluster(cl)
# stopImplicitCluster()
future:::ClusterRegistry("stop")
accuracy_rates_simu3 <- results
# accuracy_mles_simu3 <- results$accuracy_mles_simu3

# colnames(accuracy_mles_simu3) <- params_names_latex
# 
# for (param_name in c(paste0("lambda", 1:m), paste0("gamma", 1:(m ^ 2)))) {
#   sapply(accuracy_mles_simu3[accuracy_mles_simu3$par_name == param_name, 5:14], sd, na.rm = TRUE)
# }
# q <- apply(accuracy_rates_simu3,
#            2,
#            quantile.colwise)
# conf_int_simu3$Bootstrap.L[which(conf_int_simu3$m == m)] <- q[1, ]
# conf_int_simu3$Bootstrap.U[which(conf_int_simu3$m == m)] <- q[2, ]


# The profile CIs may not be sorted, so we sort them manually
# for (i in 1:length(conf_int_simu3[, 1])) {
#   row <- conf_int_simu3[i, c("Profile.L", "Profile.U")]
#   conf_int_simu3[i, c("Profile.L", "Profile.U")] <- cbind(min(row), max(row))
# }
# conf_int_simu3$m <- as.integer(conf_int_simu3$m)

# estim_benchmarks_df_simu3$m <- factor(estim_benchmarks_df_simu3$m,
#                                      levels = M_LIST_SIMU3)

# Fixes -------------------------
if (BOOTSTRAP_SAMPLES != 0) {
  # Reorder the TPM row-wise instead of column-wise
  # Rename parameters in bootstrap_simu3_results to TeX greek letters
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
  # len_par <- length(params_names_latex)
  # indices <- (length(bootstrap_simu3_results$m) + 1):(length(bootstrap_simu3_results$m) + len_par)
  bootstrap_simu3_results$all_natural_parameters$Parameter <- params_names_latex
  
  all_natural_parameters <- bootstrap_simu3_results$all_natural_parameters
  # Lexicographical parameter sort for gamma
  for (seed in unique(all_natural_parameters$seed)) {
    
    seed_indices <- all_natural_parameters$seed == seed
    truncated_table <- all_natural_parameters[seed_indices, ]
    
    # Lexicographical parameter sort for gamma (sort on the parameter name)
    new_gamma_indices_truncated_table <- order(truncated_table[gamma_indices, "Parameter"])
    # Replace rows by sorted rows
    all_natural_parameters[seed_indices, ][gamma_indices, ] <- truncated_table[gamma_indices, ][new_gamma_indices_truncated_table, ]
    
  }
  bootstrap_simu3_results$all_natural_parameters <- all_natural_parameters
}


notif(title = "Simu3 Robustness FINI",
      msg = "fin du fichier")
sink()