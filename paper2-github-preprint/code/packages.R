# type <- if (getOption("pkgType") == "both") "binary" else getOption("pkgType")
# # contriburl = contrib.url("https://cran.uib.no", "binary")
# contriburl = contrib.url("https://cran.rstudio.com", "binary")

# update.packages(ask = FALSE, type = type, contriburl = contriburl)
# # update.packages(ask = FALSE, type = type)

packages <- c("microbenchmark", "ggplot2",
              "markovchain",
              "ggthemes", "knitr", "ggtext",
              "R.utils", "lubridate", "readr", "tidyverse", "ggpubr",
              # "doParallel",
              "doFuture", "iterators", "doRNG",
              "foreach", "data.table", "itertools")
OPTIMIZER_PACKAGES <- c("TMB", "optimr", "marqLevAlg", "minqa", "BB", "purrr")
extra <- c("beepr", "ucminf", "httr", "flock")
latex <- c("xtable", "xfun", "huxtable")
packages <- c(packages, extra, OPTIMIZER_PACKAGES, latex)
if(!all( packages %in% (.packages()) )) {
  for(pkg in packages) {
    if( !require(pkg, character.only = TRUE) ) {
      # install.packages(pkg, quiet = TRUE)
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}
# "future" and "doFuture" packages get 2 identical entries in bib file -> remove one
packages <- (.packages())
packages <- packages[ packages != "doFuture" ]
knitr::write_bib(packages, 'packages.bib')

# DEBUG TMB .cpp files
# https://github.com/kaskr/TMB_contrib_R
