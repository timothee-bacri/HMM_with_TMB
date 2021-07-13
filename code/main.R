# sink("log.txt") # Redirect all messages from the console to log.txt
source("code/setup_parameters.R")
# options("warn")
# options(warn = 2)
# options(warn = 0) # DEFAULT

# knitr::knit2pdf("paper.rnw")

# Which dataset are tests run on
RUN_LAMB <- FALSE
RUN_SIMULATION1 <- FALSE
RUN_SIMULATION2 <- FALSE
# RUN_HOSPITAL <- FALSE
RUN_TINNITUS <- FALSE


begin <- Sys.time()
## ---- Lamb
if (RUN_LAMB) {
  source("code/poi_hmm_lamb.R")
  save(M_LIST_LAMB, estim_benchmarks_df_lamb, method_comparison_df_lamb, coverage_skips_lamb,
       mllk_values_lamb, mllk_times_df_lamb, conf_int_lamb,
       file = "data/results_lamb.RData")
} else {
  load("data/results_lamb.RData")
}
lambtime <- Sys.time()
# sms(paste("lamb ", round(lambtime - begin, 1), units(lambtime - begin)))
# print(lambtime - begin)
# ---- Simulation
if (RUN_SIMULATION1) {
  source("code/poi_hmm_simu1.R")
  save(M_LIST_SIMU1, estim_benchmarks_df_simu1, method_comparison_df_simu1, coverage_skips_simu1,
       mllk_values_simu1, mllk_times_df_simu1, conf_int_simu1,
       file = "data/results_simu1.RData")
} else {
  load("data/results_simu1.RData")
}
simu1time <- Sys.time()
# sms(paste("simu1 ", round(simu1time - lambtime, 1), units(simu1time - lambtime)))
if (RUN_SIMULATION2) {
  source("code/poi_hmm_simu2.R")
  save(M_LIST_SIMU2, estim_benchmarks_df_simu2, method_comparison_df_simu2, coverage_skips_simu2,
       mllk_values_simu2, mllk_times_df_simu2, conf_int_simu2,
       file = "data/results_simu2.RData")
} else {
  load("data/results_simu2.RData")
}
simu2time <- Sys.time()
# sms(paste("simu2 ", round(simu2time - simu1time, 1), units(simu2time - simu1time)))
# print(simu2time - lambtime)
## ---- Hospital
# if (RUN_HOSPITAL) {
#   source("code/poi_hmm_hosp.R")
#   save(M_LIST_HOSP, estim_benchmarks_df_hosp, method_comparison_df_hosp,
#        mllk_values_hosp, mllk_times_df_hosp, conf_int_hosp,
#        file = "data/results_hosp.RData")
# } else {
#   load("data/results_hosp.RData")
# }
# hosptime <- Sys.time()
# print(hosptime - simutime)
## ---- Tinnitus
if (RUN_TINNITUS) {
  source("code/poi_hmm_tinn.R")
  save(M_LIST_TINN, estim_benchmarks_df_tinn, method_comparison_df_tinn, coverage_skips_tinn,
       mllk_values_tinn, mllk_times_df_tinn, conf_int_tinn, consistency_estim_benchmarks_df_tinn,
       file = "data/results_tinn.RData")
} else {
  load("data/results_tinn.RData")
}
tinntime <- Sys.time()
# sms(paste("tinn ", round(tinntime - simu2time, 1), units(tinntime - simu2time)))
# print(tinntime - hosptime)
# for(a in 1:3){
#   i=0
#   message(paste0("a=",a))
#   repeat{
#     message(i)
#     i=i+1
#     if(i!=7){
#       next
#     }
#     break
#   }
#   message("out of inner loop")
# }
# sink() # Stop the log redirection
# sms(paste("lamb ", round(lambtime - begin, 1), units(lambtime - begin)))
# sms(paste("simu1 ", round(simu1time - lambtime, 1), units(simu1time - lambtime)))
# sms(paste("simu2 ", round(simu2time - simu1time, 1), units(simu2time - simu1time)))
# # sms(paste("hosp ", round(hosptime - simu2time, 1), units(hosptime - simu2time)))
# sms(paste("tinn ", round(tinntime - simu2time, 1), units(tinntime - simu2time)))
# sms(paste("total ", round(tinntime - begin, 1), units(tinntime - begin)))
# sms("FINI")
# MANUALLY CAP TMB CIs until the main code is run again
# conf_int_tinn$TMB.L <- sapply(conf_int_tinn$TMB.L, max, 0)
# conf_int_tinn$TMB.U[-1:-2] <- sapply(conf_int_tinn$TMB.U[-1:-2], min, 1)

# 
# m <- 2
# if (m == 2) {
#   true_gamma <- matrix(c(0.95, 0.05,
#                          0.15, 0.85), byrow = TRUE, nrow = m, ncol = m)
# } else if (m == 3) {
#   true_gamma <- matrix(c(0.95, 0.025, 0.025,
#                          0.05, 0.90, 0.05,
#                          0.075, 0.075, 0.85), byrow = TRUE, nrow = m, ncol = m)
# }
# true_lambda <- seq(1, 7, length.out = m)
# true_delta <- stat.dist(true_gamma)
# conf_int_simu1[conf_int_simu1$m==m, "True.value"] <- as.numeric(c(true_lambda, true_gamma, true_delta))
# m <- 3
# if (m == 2) {
#   true_gamma <- matrix(c(0.95, 0.05,
#                          0.15, 0.85), byrow = TRUE, nrow = m, ncol = m)
# } else if (m == 3) {
#   true_gamma <- matrix(c(0.95, 0.025, 0.025,
#                          0.05, 0.90, 0.05,
#                          0.075, 0.075, 0.85), byrow = TRUE, nrow = m, ncol = m)
# }
# true_lambda <- seq(1, 7, length.out = m)
# true_delta <- stat.dist(true_gamma)
# conf_int_simu2[conf_int_simu2$m==m, "True.value"] <- as.numeric(c(true_lambda, true_gamma, true_delta))


# # TEST
# conf_int_simu1[4:5, ] <- conf_int_simu1[5:4, ]
# conf_int_simu1[1:2, -2] <- conf_int_simu1[2:1, -2]
# conf_int_simu1[1:8, ]
# m=2
# lambda_indices <- 1:m
# gamma_indices <- m + 1:(m ^ 2)
# delta_indices <- m ^ 2 + m + (1:m)