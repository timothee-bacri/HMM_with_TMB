packages <- c("TMB", "markovchain", "microbenchmark", "ggplot2", "optimr",
              "ggthemes","label.switching", "readxl", "knitr", "xtable",
              "R.utils", "lubridate", "readr", "tidyverse", "ggpubr", "marqLevAlg",
              "doParallel", "foreach")
additional <- c("beepr", "devtools", "httr")
all_packages <- c(packages, additional)
if(!all( all_packages %in% (.packages()) )) {
  for(pkg in all_packages) {
    if( !require(pkg, character.only = TRUE) ) {
      install.packages(pkg, type = "binary")
      library(pkg, character.only = TRUE)
    }
  }
  
  # if(!require(notifier)) {
  #   devtools::install_version("notifier")
  #   library("notifier")
  # }
}

