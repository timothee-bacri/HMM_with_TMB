library(TMB) 
# precompile(clean = TRUE)
# precompile(PKG_LIBS = "-install_name `pwd`/$@")

compile(""norm_hmm.cpp"", "-O1 -g", DLLFLAGS = "")
gdbsource(norm_hmm.R", TRUE)
