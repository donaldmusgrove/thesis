The following R packages are required to run the included code:
gridExtra
inline
lattice
RcppArmadillo

This folder includes the following files:

1-Parcellation.R
  - Parcellates a 4D fMRI dataset into parcels with sizes that are user-specified 
  - Requires the dataset be a 4-dimensional array

2-SSGLMMonParcels.R
  - Carries out SSGLMM on each parcel in parallel

3-Plots.R
  - Extra functions for plotting results
  
################################################################################
Code folder
################################################################################
cppFunctions.cpp
  - Various C++ functions for carrying out SSGLMM, including batchmeans, and a
    truncated normal sampler
	
parcel_functions.R
  - Helper functions used to parcellate the brain

SSGLMM.R
  - Main SSGLMM function

  
################################################################################
Data folder
################################################################################
design.csv
  - Design matrix (no header) with the following columns:
    1   Fixation   
    2   Shapes   
    3   Faces
    4-9 Motion correction parameters     
    10: Ones (intercept) 

fmri4D.RData
  - Sample fMRI dataset. The first 50 time points of the fMRI dataset analyzed
    by Musgrove et al (2015) stored as a 4-dimensional array 

mask3D.RData
  - A 3-dimensional array with the binary mask of voxels to analyze

parcellatedData.RData
 - A list named parcelList with 90 entries, each of which is a parcel
   corresponding to the dataset
