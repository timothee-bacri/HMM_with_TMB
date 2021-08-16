packages <- c("TMB", "markovchain", "microbenchmark", "ggplot2", "optimr",
              "ggthemes", "knitr", "xtable",
              "R.utils", "lubridate", "readr", "tidyverse", "ggpubr", "marqLevAlg",
              "doParallel", "foreach")
if(!all( packages %in% (.packages()) )) {
  for(pkg in packages) {
    if( !require(pkg, character.only = TRUE) ) {
      install.packages(pkg, type = "binary")
      library(pkg, character.only = TRUE)
    }
  }
}