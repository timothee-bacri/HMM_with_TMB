# To produce the pdf, run the command:
# knitr::knit2pdf("paper2.rnw")

source("code/setup_parameters.R")

# Simulation
if (RUN_SIMULATION1_PARALLEL) {
  source("code/setup_parameters.R")
  lambtime_par <- Sys.time()
  source("code/poi_hmm_simu1_parallel.R", echo = TRUE, max.deparse.length = 500)
  save(M_LIST_SIMU1, benchmark_optimizer_comparison_df_simu1, optimizers_estimates_simu1, nll_simu1,
       accuracy_rates_simu1,
       bootstrap_simu1_results, simu1_data, simu1_states,
       simu1_true_lambda, simu1_true_gamma, simu1_true_delta,
       file = "data/results_simu1_parallel.RData")
  simu1time_par <- Sys.time()
  notif(msg = paste0(round(simu1time_par - lambtime_par, 1),
                     units(simu1time_par - lambtime_par),
                     "\nbenchmarks=", BENCHMARK_SAMPLES,
                     "\naccuracy_samples_amount=", length(accuracy_rates_simu1$m)),
        title = "SIMU1_PARALLEL is over")
} else {
  load("data/results_simu1_parallel.RData")
}

if (RUN_SIMULATION3_PARALLEL) {
  source("code/setup_parameters.R")
  simu1time_par <- Sys.time()
  source("code/norm_hmm_simu3_parallel.R", echo = TRUE, max.deparse.length = 500)
  save(M_LIST_SIMU3, benchmark_optimizer_comparison_df_simu3, optimizers_estimates_simu3, nll_simu3,
       accuracy_rates_simu3, bootstrap_simu3_results,
       simu3_data, simu3_states,
       simu3_true_mu, simu3_true_sigma, simu3_true_gamma, simu3_true_delta,
       simu3_parameter_possibilities_initial_amounts,
       file = "data/results_simu3_parallel.RData")
  simu3time_par <- Sys.time()
  notif(msg = paste0(round(simu3time_par - simu1time_par),
                     units(simu3time_par - simu1time_par),
                     "\nbenchmarks=", BENCHMARK_SAMPLES,
                     "\naccuracy_samples_amount=", length(accuracy_rates_simu3$m)),
        title = "SIMU3_PARALLEL is over")
} else {
  load("data/results_simu3_parallel.RData")
}

# Tinnitus
if (RUN_TINNITUS_PARALLEL) {
  source("code/setup_parameters.R")
  simu3time_par <- Sys.time()
  source("code/poi_hmm_tinn_parallel.R", echo = TRUE, max.deparse.length = 500)
  save(M_LIST_TINN, benchmark_optimizer_comparison_df_tinn, optimizers_estimates_tinn, nll_tinn,
       accuracy_rates_tinn, bootstrap_tinn_results,
       file = "data/results_tinn_parallel.RData")
  tinntime_par <- Sys.time()
  notif(msg = paste0(round(tinntime_par - simu3time_par),
                     units(tinntime_par - simu3time_par),
                     "\naccuracy_samples_amount=", length(accuracy_rates_tinn$m)),
        title = "TINN_PARALLEL is over")
} else {
  load("data/results_tinn_parallel.RData")
}