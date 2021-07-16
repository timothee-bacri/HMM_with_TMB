packages <- c("TMB", "microbenchmark", "optimr", "knitr", "xtable", "R.utils", "lubridate", "marqLevAlg", "doParallel", "foreach")

if(!all( packages %in% (.packages()) )) {
  for(pkg in all_packages) {
    if( !require(pkg, character.only = TRUE) ) {
      install.packages(pkg, type = "binary")
      library(pkg, character.only = TRUE)
    }
  }
}

