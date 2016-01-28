##Required packages

The R and C++ code contained herein include several functions for implementing the copCS model for areal data as proposed in Chapter 3 of my thesis. The C++ code makes use of the Armadillo linear algebra library, avaible via the `RcppArmadillo` R package. We also make use of existing graphical packages, including `gRim` and `RBGL`. Numerical calculations of Hessian matrices and gradient vectors are carried out using the `numDeriv` package.

```R
### First, install RBGL from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")

### Install remaining packages from CRAN
install.packages(c("gRim", "numDeriv", "RcppArmadillo"))
```

copCS.R

  - Main file: sources code used and gives an example of running the model on existing data.


##Code folder
IPS.cpp
  - Various C++ functions for carrying out IPS
	
helper_functions.R
  - Functions for fitting the copCS file


###Source an example dataset

```R
source("http://www.biostat.umn.edu/~hodges/RPLMBook/Datasets/10_Slovenian_stomach_cancer/Slovenia_stomach_cancer_data.txt")
```
