##Required packages
The R and C++ code contained herein include several functions for implementing the sparse spatial generalized linear mixed model for fMRI data as proposed in Musgrove et. al. (2015). The C++ code makes use of the Armadillo linear algebra library, avaible via the `RcppArmadillo` R package:

```R
install.packages("RcppArmadillo")
```

This folder includes the following files:

SSGLMMonParcels.R
  - Carries out SSGLMM on each parcel in parallel


##Code folder
cppFunctions.cpp
  - Various C++ functions for carrying out SSGLMM, including batchmeans, and a truncated normal sampler
	
SSGLMM.R
  - Main SSGLMM function

  
##Data folder
design.csv
- Design matrix (no header) with the following columns:
  * 1: Fixation   
  * 2: Shapes   
  * 3: Faces
  * 4-9: Motion correction parameters     
  * 10: Ones (intercept) 


parcellatedData.RData
- A list named `parcelList` with 90 entries, each of which is a parcel corresponding to the dataset


##References
```
Musgrove, D. R., Hughes, J., & Eberly, L. E. (2015). Fast, fully Bayesian spatiotemporal inference for fMRI data. Biostatistics, kxv044.
```
