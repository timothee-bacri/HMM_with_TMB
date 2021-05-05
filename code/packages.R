packages <- c("TMB", "markovchain", "microbenchmark", "ggplot2", "optimr",
              "ggthemes","label.switching", "readxl", "knitr", "xtable",
              "R.utils", "lubridate", "readr", "tidyverse", "ggpubr", "marqLevAlg")
additional <- c("beepr", "devtools", "httr")
if(!all( c(packages, additional) %in% (.packages()) )) {
  for(pkg in c(packages, additional)) {
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

