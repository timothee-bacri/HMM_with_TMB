packages <- c("TMB", "microbenchmark", "ggthemes", "knitr", "xtable",
              "R.utils", "ggpubr", "doParallel", "foreach", "mvtnorm")
if(!all( packages %in% (.packages()) )) {
  for(pkg in packages) {
    if( !require(pkg, character.only = TRUE) ) {
      install.packages(pkg, type = "binary")
      library(pkg, character.only = TRUE)
    }
  }
}
