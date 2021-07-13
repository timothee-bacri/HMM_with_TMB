rm(list = ls())

source("code/packages.R")
source("functions/utils.R")

# Sets the version R runs on and the seed it will use
# The seed used can produce different random numbers depending on the version of R
RNGversion("3.6.0")

# Used for Wald confidence intervals, approximately 1.96
q95_norm <- qnorm(1 - 0.05 / 2)

# We use nlminb because it's a fast method that works properly
# The available methods in optim package are BFGS, CG, Nelder-Mead, L-BFGS-B, nlm, nlminb,
# Rcgmin, Rvmmin, hjn
# Rvmmin and Rcgmin don't work in our case, so we don't use them

# Number of benchmark runs used to time the estimation of HMM parameters (CAN = 0)
BENCHMARK_SAMPLES <- 200

# Number of bootstrap samples used to obtain confidence intervals (CAN NOT = 1, will pose issues with parallelization)
BOOTSTRAP_SAMPLES <- 1000

# Number of confidence intervals used to obtain coverage probabilities (CAN = 0)
COVERAGE_SAMPLES <- 1000

# Number of benchmarks checking timing reliability (CAN = 0)
CONSISTENCY_BENCHMARK_TINN <- 200

# Number of hidden states used for timing HMMs on each dataset
M_LIST_LAMB <- 2
M_LIST_SIMU1 <- 2
M_LIST_SIMU2 <- 3
M_LIST_TINN <- 2

# Parameters returned by TMB_estimate
PARAMS_NAMES <- c("lambda", "gamma", "delta")

PROCEDURES <- c("DM", "TMB", "TMB_G", "TMB_H", "TMB_GH")
PROCEDURES_METHOD <- c("BFGS", "L-BFGS-B", "nlm", "nlminb", "hjn", "marqLevAlg")

# Loading the lamb dataset
load("data/fetal-lamb.RData")
lamb_data <- lamb
rm(lamb)

# Loading the tinnitus dataset, arousal for 14th patient
load("data/tinnitus.RData")

# Setting the size of the datasets
DATA_SIZE_LAMB <- length(lamb_data)
DATA_SIZE_SIMU1 <- 2000
DATA_SIZE_SIMU2 <- 5000
DATA_SIZE_TINN <- DATA_SIZE_TYT <- length(tinn_data)

# Container for the quality of fit values of HMMs for each dataset
mllk_values_lamb <- data.frame(m = integer(),
                               procedure = factor(levels = PROCEDURES),
                               mllk = numeric(),
                               AIC = numeric(),
                               BIC = numeric())
mllk_values_simu1 <- mllk_values_simu2 <- mllk_values_lamb
mllk_values_tinn <- mllk_values_lamb

# Container for the time benchmarks of HMMs parameter estimation on different optimizers for each dataset
method_comparison_df_lamb <- data.frame(time = numeric(),
                                        m = factor(levels = M_LIST_LAMB),
                                        procedure = factor(levels = PROCEDURES_METHOD),
                                        dataset_number = integer())
method_comparison_df_tinn <- data.frame(time = numeric(),
                                        m = factor(levels = M_LIST_TINN),
                                        procedure = factor(levels = PROCEDURES_METHOD),
                                        dataset_number = integer())
method_comparison_df_simu1 <- data.frame(time = numeric(),
                                         m = factor(levels = M_LIST_SIMU1),
                                         procedure = factor(levels = PROCEDURES_METHOD),
                                         dataset_number = integer())
method_comparison_df_simu2 <- data.frame(time = numeric(),
                                         m = factor(levels = M_LIST_SIMU2),
                                         procedure = factor(levels = PROCEDURES_METHOD),
                                         dataset_number = integer())

# Container for the time benchmarks of HMMs' negative likelihood function for each dataset
mllk_times_df_lamb <- data.frame(time = numeric(),
                                 m = factor(levels = M_LIST_LAMB),
                                 procedure = factor(levels = PROCEDURES))
mllk_times_df_tinn <- data.frame(time = numeric(),
                                 m = factor(levels = M_LIST_TINN),
                                 procedure = factor(levels = PROCEDURES))
mllk_times_df_simu1 <- data.frame(time = numeric(),
                                  m = factor(levels = M_LIST_SIMU1),
                                  procedure = factor(levels = PROCEDURES))
mllk_times_df_simu2 <- data.frame(time = numeric(),
                                  m = factor(levels = M_LIST_SIMU2),
                                  procedure = factor(levels = PROCEDURES))

# Container for the time benchmarks of HMMs parameter estimation for each dataset
estim_benchmarks_df_lamb <- data.frame(time = numeric(),
                                       m = integer(),
                                       procedure = factor(levels = PROCEDURES_METHOD),
                                       iterations = integer(),
                                       dataset_number = integer())
estim_benchmarks_df_simu1 <- estim_benchmarks_df_simu2 <- estim_benchmarks_df_tinn <- estim_benchmarks_df_lamb


consistency_estim_benchmarks_df_tinn <- data.frame(time = numeric(),
                                                   m = factor(levels = M_LIST_TINN),
                                                   procedure = factor(levels = PROCEDURES),
                                                   iterations = integer())

# Container for the parameters estimated by HMMs and their confidence intervals for each dataset
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
conf_int_simu1 <- conf_int_simu2 <- data.frame(m = integer(),
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
                                  "marqLevAlg_null" = 0,
                                  "marqLevAlg_converge" = 0,
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
                                   "marqLevAlg_null" = 0,
                                   "marqLevAlg_converge" = 0,
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
                                   "marqLevAlg_null" = 0,
                                   "marqLevAlg_converge" = 0,
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
                                  "marqLevAlg_null" = 0,
                                  "marqLevAlg_converge" = 0,
                                  "NA_value" = 0,
                                  "profile" = 0)


# Ensure that the algorithms BFGS and L-BFGS-B don't stop simply because of a low iteration limit
ctrl = list(maxit = 10000)

#TMB SETUP
TMB::compile("code/poi_hmm.cpp")
dyn.load(dynlib("code/poi_hmm"))
TMB::compile("code/linreg.cpp")
dyn.load(dynlib("code/linreg"))
# Optional debug feature for TMB, requires manual input immediately after being run
# TMB:::setupRStudio()