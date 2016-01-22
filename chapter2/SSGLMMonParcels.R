setwd("D:/Dropbox/Dissertation/code/ch2")


### Load example data
load("data/parcellatedData.RData")


### Load some helper functions
source("code/SSGLMM.R")
Rcpp::sourceCpp('code/cppFunctions.cpp')


### Load design matrix - Column format:
X  <- as.matrix(read.csv("data/design.csv", header=FALSE))
T  <- nrow(X)
Xc <- X[,c(4:10)]
X  <- X[,2:3]


### Example usage on a single parcel/dataset
Y  <- parcelList[[1]]$Y
A  <- parcelList[[1]]$A
tt <- main(iterations=5e2, Y, A, X, Xc, eigvecs=50)