# hmm-tmb
HMM likelihood optimization with TMB

Running all the tests:
Open the file hmm-tmb.Rproj in the main folder, or setup the working directory to the main folder manually.
Go in the folder code.
Open the file setup_parameters.R.
Change the variables RUN_LAMB, RUN_SIMULATION, and RUN_HOSPITAL to TRUE depending on what dataset you're interested in..
Run the file main.R.

Accessing the data:
Open the file hmm-tmb.Rproj in the main folder, or setup the working directory to the main folder manually.
Go in the folder code.
Run the file main.R.

Explanation of the files:
_Folder code
main.R can load all the data or can run the code and save the results.
packages.R installs all necessary packages, make sure it runs without error.
linreg.cpp is the C++ code specifying the objective function for the linear regression model.
linreg.dll and linreg.o are produced by linreg.cpp and are used by TMB.
poi_hmm.cpp poi_hmm.dll poi_hmm.o are similar but used for the Poisson HMM.
poi_hmm_*****.R run the timing procedures for the 3 datasets.
_Folder data
data_*****.RData contain the results of the tests. See code/main.R and code/setup_parameters.R for more details.
grouped_hospital_data.csv contains the hospital dataset.
fetal-lamb.RData and leroux.txt contain the lamb dataset. We use the first one for convenience reasons.
_Folder functions
utils.cpp contains functions used with TMB inside the C++ code code/poi_hmm.cpp.
utils.R contains functions used in R to process the parameters and the results.
utils_linreg.cpp contains functions used with TMB inside the C++ code code/linreg.cpp.