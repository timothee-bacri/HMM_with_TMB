# sink("log.txt") # Redirect all messages from the console to log.txt
source("code/setup_parameters.R")

# knitr::knit2pdf("paper.rnw")

# Which dataset are tests run on
RUN_LAMB <- FALSE
RUN_SIMULATION1 <- FALSE
RUN_SIMULATION2 <- FALSE
RUN_TINNITUS <- FALSE

## ---- Lamb
if (RUN_LAMB) {
  source("code/poi_hmm_lamb.R")
  save(M_LIST_LAMB, estim_benchmarks_df_lamb, method_comparison_df_lamb, coverage_skips_lamb,
       mllk_values_lamb, mllk_times_df_lamb, conf_int_lamb,
       file = "data/results_lamb.RData")
} else {
  load("data/results_lamb.RData")
}

# ---- Simulation
if (RUN_SIMULATION1) {
  source("code/poi_hmm_simu1.R")
  save(M_LIST_SIMU1, estim_benchmarks_df_simu1, method_comparison_df_simu1, coverage_skips_simu1,
       mllk_values_simu1, mllk_times_df_simu1, conf_int_simu1,
       file = "data/results_simu1.RData")
} else {
  load("data/results_simu1.RData")
}

if (RUN_SIMULATION2) {
  source("code/poi_hmm_simu2.R")
  save(M_LIST_SIMU2, estim_benchmarks_df_simu2, method_comparison_df_simu2, coverage_skips_simu2,
       mllk_values_simu2, mllk_times_df_simu2, conf_int_simu2,
       file = "data/results_simu2.RData")
} else {
  load("data/results_simu2.RData")
}
simu2time <- Sys.time()

## ---- Tinnitus
if (RUN_TINNITUS) {
  source("code/poi_hmm_tinn.R")
  save(M_LIST_TINN, estim_benchmarks_df_tinn, method_comparison_df_tinn, coverage_skips_tinn,
       mllk_values_tinn, mllk_times_df_tinn, conf_int_tinn, consistency_estim_benchmarks_df_tinn,
       file = "data/results_tinn.RData")
} else {
  load("data/results_tinn.RData")
}
