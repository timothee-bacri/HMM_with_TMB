sink("log.txt") # Redirect all messages from the console to log.txt
source("code/setup_parameters.R")
begin <- Sys.time()
## ---- Lamb
if (RUN_LAMB) {
  source("code/poi_hmm_lamb.R")
  save(M_LIST_LAMB, estim_benchmarks_df_lamb, method_comparison_df_lamb,
       mllk_values_lamb, mllk_times_df_lamb, conf_int_lamb, consistency_estim_benchmarks_df_lamb,
       file = "data/results_lamb.RData")
} else {
  load("data/results_lamb.RData")
}
lambtime <- Sys.time()
# print(lambtime - begin)
## ---- Simulation
if (RUN_SIMULATION) {
  source("code/poi_hmm_simu.R")
  save(M_LIST_SIMU, estim_benchmarks_df_simu, method_comparison_df_simu,
       mllk_values_simu, mllk_times_df_simu, conf_int_simu, consistency_estim_benchmarks_df_simu,
       file = "data/results_simu.RData")
} else {
  load("data/results_simu.RData")
}
simutime <- Sys.time()
# print(simutime - lambtime)
## ---- Hospital
if (RUN_HOSPITAL) {
  source("code/poi_hmm_hosp.R")
  save(M_LIST_HOSP, estim_benchmarks_df_hosp, method_comparison_df_hosp,
       mllk_values_hosp, mllk_times_df_hosp, conf_int_hosp, consistency_estim_benchmarks_df_hosp,
       file = "data/results_hosp.RData")
} else {
  load("data/results_hosp.RData")
}
hosptime <- Sys.time()
# print(hosptime - simutime)
## ---- Tinnitus
if (RUN_TINNITUS) {
  source("code/poi_hmm_tinn.R")
  save(M_LIST_TINN, estim_benchmarks_df_tinn, method_comparison_df_tinn,
       mllk_values_tinn, mllk_times_df_tinn, conf_int_tinn, consistency_estim_benchmarks_df_tinn,
       file = "data/results_tinn.RData")
} else {
  load("data/results_tinn.RData")
}
tinntime <- Sys.time()
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
sink() # Stop the log redirection
# sms(paste("lamb ", round(lambtime - begin, 1), units(lambtime - begin)))
# sms(paste("simu ", round(simutime - lambtime, 1), units(simutime - lambtime)))
# sms(paste("hosp ", round(hosptime - simutime, 1), units(hosptime - simutime)))
# sms(paste("tinn ", round(tinntime - hosptime, 1), units(tinntime - hosptime)))
# sms(paste("total ", round(tinntime - begin, 1), units(tinntime - begin)))
