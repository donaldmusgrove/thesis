################################################################################
### Main function call
################################################################################
main <- function(iterations, Y, A, X, Xc, eigvecs){

  ### Set up constants
  np   <- ncol(Y)
  p    <- ncol(X)
  tv   <- length(Y[,1])
  u    <- ncol(Xc)
  w    <- diag(tv) - Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)
  www  <- eigen(w)$vector[,1:(tv-u)]
  Y    <- apply(Y,2,function(x) t(www)%*%x)
  X    <- t(www)%*%X
  tv   <- nrow(Y)

  ### Setup spatial dependency structure
  A    <- A
  Q    <- diag(c(A%*%rep(1,np)),nrow=np)-A
  vvv  <- eigen( A )
  v1   <- range(vvv$values)
  v2   <- 2/(v1[2]-v1[1])
  www  <- -1 + (vvv$values-v1[1] )*v2
  M    <- vvv$vector[,1:sum(www>0.05)]
  if(eigvecs != 0){
    M   <- vvv$vectors[,1:eigvecs]
  }
  Qs  <- t(M)%*%Q%*%M
  q   <- ncol(M)
  Sigma     <- chol2inv(chol(Qs+t(M)%*%M))
  Sigmachol <- t(chol(Sigma))
  nu        <- diag( diag(np) + M%*%solve(Qs)%*%t(M) )

  ###Initialize the chain(s)
  beta0              <- t(sapply(1:np, function(i) summary(lm(Y[,i] ~ X))$coef[,c(1,4)] ))
    beta.chain       <- array(,dim=c(np,p,1))
    beta.chain[,,1]  <- beta0[,1:p]
  gamma.chain        <- array(,dim=c(np,p,iterations))
    gamma.chain[,,1] <- ifelse(beta0[,(p+1):(2*p)] < 0.05,1,0)
    beta.chain[,,1]  <- sapply(1:p, function(j) beta.chain[,j,1]*gamma.chain[,j,1])
  sigma2.chain       <- matrix(0,np,1)
    sigma2.chain[,1] <- sapply(1:np, function(i) sum((Y[,i] - X%*%beta.chain[i,,1])^2)/tv)
  rho.chain          <- array(,dim=c(np,2,1))
    for(i in 1:np){
      w   <- Y[,i] - X%*%beta.chain[i,,1]
      wL  <- cbind(w[2:(tv-1)], w[1:(tv-2)])
      rho.chain[i,,1] <- solve(t(wL)%*%wL)%*%t(wL)%*%w[3:tv]
    }
  tau2.chain         <- matrix(0,p,1)
    tau2.chain[,1]   <- sapply(1:p, function(j) 1/(sum(gamma.chain[,j,1]))*sum(gamma.chain[,j,1]*beta.chain[,j,1]^2) )
  phi.chain          <- array(,dim=c(q,p,1))
    phi.chain[,,1]   <- sapply(1:p, function(i) glm(gamma.chain[,i,1] ~ M-1)$coef)
  theta.chain        <- array(,dim=c(np,p,1))
    theta.chain[,,1] <- matrix(pnorm(qnorm(0.02) + M%*%phi.chain[,,1]),ncol=p)
  kappa.chain        <- matrix(0,p,1)
    kappa.chain[,1]  <- sapply(1:p, function(i) q/t(phi.chain[,i,1])%*%Qs%*%phi.chain[,i,1] )
 
  ### MCMC sampling
  for(k in 2:iterations){
    eta               <- pnorm(qnorm(0.02) + theta.chain[,,1])
    reg.out           <- gammaSampler(Y, as.matrix(X), as.matrix(eta), as.matrix(beta.chain[,,1]), 
                                      as.matrix(gamma.chain[,,k-1]), rho.chain[,,1], tau2.chain[,1], sigma2.chain[,1])
    gamma.chain[,,k]  <- reg.out$gamma
    beta.chain[,,1]   <- reg.out$beta
    tau2.chain[,1]    <- tau2Sampler(X, as.matrix(beta.chain[,,1]), as.matrix(gamma.chain[,,k]), tau2.chain[,1])
    sigma2.chain[,1]  <- sigma2Sampler(Y, X, as.matrix(beta.chain[,,1]), rho.chain[,,1])
    rho.chain[,,1]    <- rhoSampler(Y, X, as.matrix(beta.chain[,,1]), rho.chain[,,1], sigma2.chain[,1])

    spat.out          <- phiSampler(as.matrix(phi.chain[,,1]), as.matrix(theta.chain[,,1]), 
                                    as.matrix(gamma.chain[,,k]), kappa.chain[,1], M, Sigma, Sigmachol, nu)
    phi.chain[,,1]    <- spat.out$phi
    theta.chain[,,1]  <- spat.out$theta
    kappa.chain[,1]   <- kappaSampler(as.matrix(theta.chain[,,1]), nu)
  }

  gammaEST        <- lapply(1:p, function(j) bmmatcpp(t(gamma.chain[,j,1:iterations])) )
  gammaEST_burnin <- lapply(1:p, function(j) bmmatcpp(t(gamma.chain[,j,(iterations/2):iterations])) )
  rm(gamma.chain)

  res             <- list(gammaEST=gammaEST, gammaEST_burnin=gammaEST_burnin)

  res
}