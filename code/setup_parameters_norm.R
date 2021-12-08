rm(list = ls())

source("./code/packages.R")
source("./functions/utils.R")
source("./functions/utils_personal.R")

# Sets the version R runs on and the seed it will use
# The seed used can produce different random numbers depending on the version of R
RNGversion("3.6.0")
set.seed(123)

# Which dataset are tests run on
RUN_SP500 <- FALSE

# Used for Wald confidence intervals, approximately 1.96
q95_norm <- qnorm(1 - 0.05 / 2)

# We use nlminb because it's a fast method that works properly
# The possible methods are BFGS, CG, Nelder-Mead, L-BFGS-B, nlm, nlminb,
# Rcgmin, Rvmmin, hjn
# Rvmmin and Rcgmin don't work in our case, so we don't use them

# Number of benchmarks used to time estimation of HMM parameters
BENCHMARK_TRIALS <- 5

# Number of bootstrap samples used to retrieve confidence intervals
BOOTSTRAP_SAMPLES <- 500

# Number of hidden states used for timing HMMs on each dataset
# M_LIST_SP500 <- c(2, 3, 4)
M_LIST_SP500 <- c(2, 3, 4)

# Parameters returned by TMB.estimate
params_names <- c("mu", "sigma", "gamma", "delta")

# Name if TMB provides gradient/hessian
#
#                        Gradient not provided   Gradient provided
# Hessian not provided          TMB1                   TMB3
# Hessian provided              TMB2                   TMB4
PROCEDURES <- c("DM", "TMB1", "TMB2", "TMB3", "TMB4")

# Loading the sp500 dataset
# SP500 <- read_csv("data/S&P500.csv") %>%
#   mutate(lag = lag(`Adj Close`),
#          lr = log(`Adj Close`/lag)) %>%
#   select(lr) %>%
#   drop_na()
# SP500 <- 100 * SP500$lr
# save(SP500, file = "./data/SP500.RData")
load("./data/SP500.RData")

# Setting the size of the datasets
DATA_SIZE_SP500 <- length(SP500)

# Container for the quality of fit values of HMMs for each dataset
mllk_values_SP500 <- data.frame(m = integer(),
                                procedure = factor(levels = PROCEDURES),
                                mllk = numeric(),
                                AIC = numeric(),
                                BIC = numeric())

# Container for the time benchmarks of HMMs parameter estimation on different optimizers for each dataset
method_comparison_SP500 = list()

# Container for the time benchmarks of HMMs' negative likelihood function for each dataset
mllk_times_SP500 <- list()

# Container for the time benchmarks of HMMs parameter estimation for each dataset
benchmarks_SP500 <- list()

# More practical container for the time benchmarks of HMMs parameter estimation for each dataset
benchmarks_df_SP500 <- data.frame(time = numeric(),
                                  m = factor(levels = M_LIST_SP500),
                                  procedure = factor(levels = PROCEDURES))

# Container for the parameters estimated by HMMs and their confidence intervals for each dataset
conf_int_SP500 <- data.frame(m = integer(),
                             Parameter = character(),
                             Parameter.estimate = numeric(),
                             Profile.L = numeric(),
                             Profile.U = numeric(),
                             Bootstrap.L = numeric(),
                             Bootstrap.U = numeric(),
                             TMB.L = numeric(),
                             TMB.U = numeric(),
                             stringsAsFactors = FALSE)

# Ensure that the algorithms BFGS and L-BFGS-B don't stop simply because of a low iteration limit
ctrl = list(maxit = 10000)

#TMB SETUP
TMB::compile("code/norm_hmm.cpp")
dyn.load(dynlib("code/norm_hmm"))
# dyn.unload("code/norm_hmm")
# Optional debug feature for TMB, requires manual input immediately after being run
# TMB:::setupRStudio()