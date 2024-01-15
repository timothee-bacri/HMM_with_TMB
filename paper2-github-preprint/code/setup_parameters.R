rm(list = ls())

# Which data sets are analyses run on
RUN_SIMULATION1_PARALLEL <- FALSE
RUN_SIMULATION3_PARALLEL <- FALSE
RUN_TINNITUS_PARALLEL <- FALSE

# Set COMPILE_AND_LOAD_TMB_FILES to TRUE for compiling the rnw file
COMPILE_AND_LOAD_TMB_FILES <- FALSE
# COMPILE_AND_LOAD_TMB_FILES <- TRUE

if (COMPILE_AND_LOAD_TMB_FILES) {
  source("code/packages.R")
}
source("functions/utils.R")
source("functions/norm_utils.R")

# Some commands contain private info and are not on GitHub
notif <- function(msg, title, priority) {}

# How many cores to use in parallel
# CORES <- parallel::detectCores() / 2
# CORES <- future::availableCores() / 2
CORES <- 4

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

# Number of hidden states used for timing HMMs on each dataset
M_LIST_SIMU1 <- 2
M_LIST_SIMU3 <- 2
M_LIST_TINN <- 2

# Percentage
NLL_THRESHOLD <- 0.05

# Parameters returned by TMB.estimate
PARAMS_NAMES <- c("lambda", "gamma", "delta")
PARAMS_NAMES_NORM <- c("mu", "sigma", "gamma", "delta")

OPTIMIZERS_METHOD <- c("BFGS", "BFGS_gr",
                       "L-BFGS-B", "L-BFGS-B_gr",
                       "CG", "CG_gr",
                       "Nelder-Mead",
                       "nlm", "nlm_gr", "nlm_he", "nlm_grhe",
                       "nlminb", "nlminb_gr", "nlminb_he", "nlminb_grhe",
                       "hjn",
                       "marqLevAlg", "marqLevAlg_gr", "marqLevAlg_he", "marqLevAlg_grhe",
                       "newuoa",
                       "BBoptim")
# "BBoptim_gr")
OPTIMIZERS_METHOD_UNDERSCORE <- gsub(pattern = "-", replacement = "_", x = OPTIMIZERS_METHOD)
OPTIMIZERS_METHOD_BACKSLASH <- gsub(pattern = "_", replacement = "\\\\_", x = OPTIMIZERS_METHOD)
OPTIMIZERS_MISSING_ITERATIONS <- c("BFGS", "L-BFGS-B", "CG", "Nelder-Mead", "hjn", "newuoa")
OPTIMIZERS_MISSING_ITERATIONS_FULL <- c(OPTIMIZERS_MISSING_ITERATIONS, "BFGS_gr", "L-BFGS-B_gr", "CG_gr")

# Ensure that the optimizers don't stop simply because of a low iteration/function evaluation limit
max_iterations <- 10000
CONTROL_ARGS <- list(BFGS = list(trace = 0,
                                 maxit = max_iterations),
                     
                     "L-BFGS-B" = list(trace = 0,
                                       maxit = max_iterations),
                     
                     CG = list(trace = 0,
                               maxit = max_iterations),
                     
                     "Nelder-Mead" = list(trace = 0,
                                          maxit = max_iterations),
                     
                     nlm = list(iterlim = max_iterations),
                     
                     nlminb = list(iter.max = max_iterations),
                     
                     marqLevAlg = list(maxiter = max_iterations),
                     
                     BBoptim = list(trace = FALSE,
                                    maxit = max_iterations))

# Loading the tinnitus data set
load("data/tinnitus.RData")

# # Download the S&P500 data
# library(quantmod)
# SP500 <- getSymbols(Symbols = "^GSPC",
#                     src = "yahoo",
#                     auto.assign = FALSE,
#                     from = as.Date("1980-01-01"),
#                     to = as.Date("2022-09-30"))
# # # Retrieve weekly returns
# sp500_data <- periodReturn(x = SP500,
#                            type = "log",
#                            period = "weekly")
# # # First return does not exist
# sp500_data <- sp500_data[-1, ]
# sp500_dates <- as.Date(rownames(data.frame(sp500_data)))
# sp500_data <- as.vector(sp500_data)
# # write(sp500_data, "data/SP500.txt", sep = "\n")
sp500_data <- read.table("data/SP500.txt", quote="\"", comment.char="")
sp500_data <- as.numeric(unlist(sp500_data))


# Loading the hospital data set
# The DATE column requires a line of code to be treated as a proper datetime item (POSIXct)
# To that end, we use the lubridate package
# library(dplyr)
# library(lubridate)
# full_hosp_data <- as_tibble(read.csv("data/grouped_hospital_data.csv",
#                                      stringsAsFactors = FALSE)) %>%
#   mutate(DATE = ymd_hms(DATE))
# hosp_data <- full_hosp_data$PATIENTS
# 
# save(full_hosp_data, file = "data/full_hosp_data.RData")
# save(hosp_data, file = "data/hosp_data.RData")
load("data/full_hosp_data.RData")
load("data/hosp_data.RData")

# Load the soap data set
soap_data <- read.table("data/soap.txt", quote="\"", comment.char="")
soap_data <- as.numeric(unlist(soap_data))

# Setting the size of the datasets
DATA_SIZE_SIMU1 <- 200
DATA_SIZE_SIMU3 <- 200
DATA_SIZE_TINN <- DATA_SIZE_TYT <- length(tinn_data)
DATA_SIZE_HOSP <- length(hosp_data)
DATA_SIZE_SOAP <- length(soap_data)
DATA_SIZE_SP500 <- length(sp500_data)

# Container for the time benchmarks of HMMs parameter estimation on different optimizers for each dataset
benchmark_optimizer_comparison_df_tinn <- data.frame(time = numeric(),
                                                     m = factor(levels = M_LIST_TINN),
                                                     optimizer = factor(levels = OPTIMIZERS_METHOD),
                                                     dataset_number = integer())
benchmark_optimizer_comparison_df_simu1 <- data.frame(time = numeric(),
                                                      m = factor(levels = M_LIST_SIMU1),
                                                      optimizer = factor(levels = OPTIMIZERS_METHOD),
                                                      dataset_number = integer())
benchmark_optimizer_comparison_df_simu3 <- data.frame(time = numeric(),
                                                      m = factor(levels = M_LIST_SIMU3),
                                                      optimizer = factor(levels = OPTIMIZERS_METHOD),
                                                      dataset_number = integer())

# Containers for the time benchmarks of HMMs parameter estimation for each dataset
estim_benchmarks_df_simu1 <- data.frame(time = numeric(),
                                        m = integer(),
                                        optimizer = factor(levels = OPTIMIZERS_METHOD),
                                        iterations = integer(),
                                        dataset_number = integer())
estim_benchmarks_df_simu3 <- estim_benchmarks_df_simu1
estim_benchmarks_df_tinn <- estim_benchmarks_df_simu1

accuracy_mles_simu1 <- accuracy_mles_simu3 <- accuracy_mles_tinn <- data.frame()

# Robustness settings Simu1
true_lambda_1_simu1 <- seq(1, 4, length.out = M_LIST_SIMU1)
true_lambda_2_simu1 <- seq(5, 7, length.out = M_LIST_SIMU1)
true_lambda_simu1 <- list(true_lambda_1_simu1, true_lambda_2_simu1)
true_gamma_1_simu1 <- matrix(c(0.9, 0.1,
                               0.2, 0.8),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU1,
                             ncol = M_LIST_SIMU1)
true_gamma_2_simu1 <- matrix(c(0.7, 0.3,
                               0.8, 0.2),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU1,
                             ncol = M_LIST_SIMU1)
true_gamma_3_simu1 <- matrix(c(0.1, 0.9,
                               0.8, 0.2),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU1,
                             ncol = M_LIST_SIMU1)
true_gamma_4_simu1 <- matrix(c(0.55, 0.45,
                               0.45, 0.55),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU1,
                             ncol = M_LIST_SIMU1)
true_gamma_simu1 <- list(true_gamma_1_simu1,
                         true_gamma_2_simu1,
                         true_gamma_3_simu1,
                         true_gamma_4_simu1)
if (!require(purrr)) {
  install.packages("purrr")
  library(purrr)
}
true_params_simu1 <- cross2(true_lambda_simu1, true_gamma_simu1)
# Add the stationary distributions
true_params_simu1 <- lapply(true_params_simu1, function(x) {
  names(x) <- c("lambda", "gamma")
  x[["delta"]] <- stat.dist(x[[2]])
  return(x)
})


# Robustness settings Simu3
true_mu_1_simu3 <- seq(-2, 2, length.out = M_LIST_SIMU3)
true_mu_2_simu3 <- seq(-1, 4, length.out = M_LIST_SIMU3)
true_mu_simu3 <- list(true_mu_1_simu3, true_mu_2_simu3)
true_sigma_1_simu3 <- seq(1.5, 2.5, length.out = M_LIST_SIMU3)
true_sigma_2_simu3 <- seq(0.5, 1.5, length.out = M_LIST_SIMU3)
true_sigma_simu3 <- list(true_sigma_1_simu3, true_sigma_2_simu3)
true_gamma_1_simu3 <- matrix(c(0.9, 0.1,
                               0.2, 0.8),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU3,
                             ncol = M_LIST_SIMU3)
true_gamma_2_simu3 <- matrix(c(0.7, 0.3,
                               0.8, 0.2),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU3,
                             ncol = M_LIST_SIMU3)
true_gamma_3_simu3 <- matrix(c(0.1, 0.9,
                               0.8, 0.2),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU3,
                             ncol = M_LIST_SIMU3)
true_gamma_4_simu3 <- matrix(c(0.55, 0.45,
                               0.45, 0.55),
                             byrow = TRUE,
                             nrow = M_LIST_SIMU3,
                             ncol = M_LIST_SIMU3)
true_gamma_simu3 <- list(true_gamma_1_simu3,
                         true_gamma_2_simu3,
                         true_gamma_3_simu3,
                         true_gamma_4_simu3)
true_params_simu3 <- cross3(true_mu_simu3, true_sigma_simu3, true_gamma_simu3)
# Add the stationary distributions
true_params_simu3 <- lapply(true_params_simu3, function(x) {
  names(x) <- c("mu", "sigma", "gamma")
  x[["delta"]] <- stat.dist(x[[3]])
  return(x)
})


if (COMPILE_AND_LOAD_TMB_FILES) {
  # TMB SETUP
  if ("poi_hmm" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib("code/poi_hmm"))
  }
  if ("poi_hmm_smoothing.cpp" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib("code/poi_hmm_smoothing"))
  }
  if ("poi_hmm_MSE_estimators" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib("code/poi_hmm_MSE_estimators"))
  }
  if ("norm_hmm" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib("code/norm_hmm"))
  }
  if ("norm_hmm_smoothing.cpp" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib("code/norm_hmm_smoothing"))
  }
  # TMB::precompile()
  TMB::compile("code/poi_hmm.cpp")
  dyn.load(dynlib("code/poi_hmm"))
  TMB::compile("code/norm_hmm.cpp")
  dyn.load(dynlib("code/norm_hmm"))
  TMB::compile("code/poi_hmm_smoothing.cpp")
  dyn.load(dynlib("code/poi_hmm_smoothing"))
  TMB::compile("code/norm_hmm.cpp")
  dyn.load(dynlib("code/norm_hmm"))
  TMB::compile("code/norm_hmm_smoothing.cpp")
  dyn.load(dynlib("code/norm_hmm_smoothing"))
  TMB::compile("code/poi_hmm_MSE_estimators.cpp")
  dyn.load(dynlib("code/poi_hmm_MSE_estimators"))
  # Optional debug feature for TMB, requires manual input immediately after being run
  # TMB:::setupRStudio()
}

# Graphical settings for nice plots
ERRORBAR_LINE_WIDTH <- 0.2
ERRORBAR_END_WIDTH <- 1.1
POINT_SIZE <- 0.2
LINE_SIZE <- 0.25
COLORS <- c("dodgerblue2", "tomato3", "green", "purple", "gold")
