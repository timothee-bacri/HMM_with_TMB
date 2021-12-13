rm(list = ls())

source("code/packages.R")
source("functions/utils.R")

CORES <- 32

# Sets the version R runs on and the seed it will use
# The seed used can produce different random numbers depending on the version of R
RNGversion("3.6.0")

# Used for Wald confidence intervals, approximately 1.96
q95_norm <- qnorm(1 - 0.05 / 2)

# Number of benchmark runs used to time the estimation of HMM parameters (CAN = 0)
BENCHMARK_SAMPLES <- 200
# BENCHMARK_SAMPLES <- 2

# Number of bootstrap samples used to obtain confidence intervals (CAN = 0)
BOOTSTRAP_SAMPLES <- 1000
# BOOTSTRAP_SAMPLES <- 2

# Number of confidence intervals used to obtain coverage probabilities (CAN = 0)
COVERAGE_SAMPLES <- 1000
# COVERAGE_SAMPLES <- 2

# Number of benchmarks checking timing reliability (CAN = 0)
CONSISTENCY_BENCHMARK_TINN <- 200
# CONSISTENCY_BENCHMARK_TINN <- 2

# Number of hidden states used for timing HMMs on each dataset
M_LIST_LAMB <- 2
M_LIST_SIMU1 <- 2
M_LIST_SIMU2 <- 3
M_LIST_SIMU3 <- 4
M_LIST_SIMU4 <- 4
M_LIST_TINN <- 2

# Parameters returned by TMB.estimate
PARAMS_NAMES <- c("lambda", "gamma", "delta")

PROCEDURES <- c("DM", "TMB", "TMB_G", "TMB_H", "TMB_GH")

# Loading the lamb dataset
load("data/fetal-lamb.RData")
lamb_data <- lamb
rm(lamb)

# Loading the tinnitus dataset
load("data/tinnitus.RData")

# Setting the size of the datasets
DATA_SIZE_LAMB <- length(lamb_data)
DATA_SIZE_SIMU1 <- 2000
DATA_SIZE_SIMU2 <- 5000
DATA_SIZE_SIMU3 <- 2000
DATA_SIZE_SIMU4 <- 5000
DATA_SIZE_TINN <- DATA_SIZE_TYT <- length(tinn_data)

# Containers for the time benchmarks of HMMs parameter estimation for each dataset
estim_benchmarks_df_lamb <- data.frame(time = numeric(),
                                       m = integer(),
                                       procedure = factor(levels = PROCEDURES),
                                       iterations = integer(),
                                       dataset_number = integer())
estim_benchmarks_df_simu1 <- estim_benchmarks_df_simu2 <- estim_benchmarks_df_lamb
estim_benchmarks_df_simu3 <- estim_benchmarks_df_simu4 <- estim_benchmarks_df_lamb
estim_benchmarks_df_tinn <- estim_benchmarks_df_lamb

consistency_estim_benchmarks_df_tinn <- data.frame(time = numeric(),
                                                   m = factor(levels = M_LIST_TINN),
                                                   procedure = factor(levels = PROCEDURES),
                                                   iterations = integer())

# Containers for the parameters estimated by HMMs and their confidence intervals for each dataset
conf_int_lamb <- conf_int_tinn <- data.frame(m = integer(),
                                             Parameter = character(),
                                             Estimate = numeric(),
                                             TMB.L = numeric(),
                                             TMB.U = numeric(),
                                             Profile.L = numeric(),
                                             Profile.U = numeric(),
                                             Bootstrap.L = numeric(),
                                             Bootstrap.U = numeric(),
                                             Coverage.TMB = numeric(),
                                             Coverage.Profile = numeric(),
                                             Coverage.Bootstrap = numeric(),
                                             stringsAsFactors = FALSE)

# For the simulation, we know the true parameter value being estimated
conf_int_simu1 <- conf_int_simu2 <- conf_int_simu3 <- conf_int_simu4 <- data.frame(m = integer(),
                                                                                   Parameter = character(),
                                                                                   True.value = numeric(),
                                                                                   Estimate = numeric(),
                                                                                   TMB.L = numeric(),
                                                                                   TMB.U = numeric(),
                                                                                   Profile.L = numeric(),
                                                                                   Profile.U = numeric(),
                                                                                   Bootstrap.L = numeric(),
                                                                                   Bootstrap.U = numeric(),
                                                                                   Coverage.TMB = numeric(),
                                                                                   Coverage.Profile = numeric(),
                                                                                   Coverage.Bootstrap = numeric(),
                                                                                   stringsAsFactors = FALSE)

# Number of datasets discarded due to an issue
coverage_skips_lamb <- data.frame("m" = M_LIST_LAMB,
                                  "state_number" = 0,
                                  "TMB_null" = 0,
                                  "TMB_converge" = 0,
                                  "TMB_G_null" = 0,
                                  "TMB_G_converge" = 0,
                                  "TMB_H_null" = 0,
                                  "TMB_H_converge" = 0,
                                  "TMG_GH_null" = 0,
                                  "TMG_GH_converge" = 0,
                                  "NA_value" = 0,
                                  "profile" = 0)
coverage_skips_simu1 <- data.frame("m" = M_LIST_SIMU1,
                                   "state_number" = 0,
                                   "TMB_null" = 0,
                                   "TMB_converge" = 0,
                                   "TMB_G_null" = 0,
                                   "TMB_G_converge" = 0,
                                   "TMB_H_null" = 0,
                                   "TMB_H_converge" = 0,
                                   "TMG_GH_null" = 0,
                                   "TMG_GH_converge" = 0,
                                   "NA_value" = 0,
                                   "profile" = 0)
coverage_skips_simu2 <- data.frame("m" = M_LIST_SIMU2,
                                   "state_number" = 0,
                                   "TMB_null" = 0,
                                   "TMB_converge" = 0,
                                   "TMB_G_null" = 0,
                                   "TMB_G_converge" = 0,
                                   "TMB_H_null" = 0,
                                   "TMB_H_converge" = 0,
                                   "TMG_GH_null" = 0,
                                   "TMG_GH_converge" = 0,
                                   "NA_value" = 0,
                                   "profile" = 0)
coverage_skips_simu3 <- data.frame("m" = M_LIST_SIMU3,
                                   "state_number" = 0,
                                   "TMB_null" = 0,
                                   "TMB_converge" = 0,
                                   "TMB_G_null" = 0,
                                   "TMB_G_converge" = 0,
                                   "TMB_H_null" = 0,
                                   "TMB_H_converge" = 0,
                                   "TMG_GH_null" = 0,
                                   "TMG_GH_converge" = 0,
                                   "NA_value" = 0,
                                   "profile" = 0)
coverage_skips_simu4 <- data.frame("m" = M_LIST_SIMU4,
                                   "state_number" = 0,
                                   "TMB_null" = 0,
                                   "TMB_converge" = 0,
                                   "TMB_G_null" = 0,
                                   "TMB_G_converge" = 0,
                                   "TMB_H_null" = 0,
                                   "TMB_H_converge" = 0,
                                   "TMG_GH_null" = 0,
                                   "TMG_GH_converge" = 0,
                                   "NA_value" = 0,
                                   "profile" = 0)
coverage_skips_tinn <- data.frame("m" = M_LIST_TINN,
                                  "state_number" = 0,
                                  "TMB_null" = 0,
                                  "TMB_converge" = 0,
                                  "TMB_G_null" = 0,
                                  "TMB_G_converge" = 0,
                                  "TMB_H_null" = 0,
                                  "TMB_H_converge" = 0,
                                  "TMG_GH_null" = 0,
                                  "TMG_GH_converge" = 0,
                                  "NA_value" = 0,
                                  "profile" = 0)


# TMB SETUP
# TMB::precompile()
TMB::compile("code/poi_hmm.cpp")
dyn.load(dynlib("code/poi_hmm"))
TMB::compile("code/norm_hmm.cpp")
dyn.load(dynlib("code/norm_hmm"))
TMB::compile("code/mvnorm_hmm.cpp")
dyn.load(dynlib("code/mvnorm_hmm"))
TMB::compile("code/linreg.cpp")
dyn.load(dynlib("code/linreg"))
# Optional debug feature for TMB, requires manual input immediately after being run
# TMB:::setupRStudio()