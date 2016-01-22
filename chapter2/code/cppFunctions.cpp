#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
// Batchmeans on a single vector
vec bmcpp(vec X){

  int n  = X.n_elem;
  int b  = floor(pow(n,0.5));
  int a  = floor(n/b);
  double xl, xu, xmean, muhat, varhat, e1;

  vec y(a);
  vec v(a);
  vec est(2);

  for (int i=1; i<(a+1); i++){
    xl      =  (i - 1) * b ;
    xu      =  i * b - 1;
    xmean   =  mean(X.rows(xl,xu));     
    y(i-1)  =  xmean;
  }

  muhat = mean(y);

  v.fill(muhat);
  v  = y - v;
  e1 = sum(trans(v) * v);

  varhat = b * e1 / (a - 1);

  est(0) = muhat;
  est(1) = pow( varhat/n, 0.5 );

  return(est) ;
}


// [[Rcpp::export]]
// Batchmeans on a matrix
mat bmmatcpp(mat M){

  int num = M.n_cols;
  mat bmvals(num,2);
  vec bmout(2);

  for (int i=0; i<num; i++){
    bmout         = bmcpp(M.col(i));
    bmvals.row(i) = trans(bmout);
  }

  return(bmvals);
}


// ======================================================================
// Truncate normal sampling functions
// ======================================================================
// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================
double norm_rs(double a, double b){
  double  x;
  x = Rf_rnorm(0.0, 1.0);
  while( (x < a) || (x > b) ) x = norm_rand();
  return x;
}


// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================
double half_norm_rs(double a, double b){
  double   x;
  x = fabs(norm_rand());
  while( (x<a) || (x>b) ) x = fabs(norm_rand());
  return x;
}


// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling. 
// ======================================================================
double unif_rs(double a, double b){
  double xstar, logphixstar, x, logu;

  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0) xstar = 0.0;
  else xstar = a;
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar)){
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}


// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================
double exp_rs(double a, double b){
  double  z, u, rate;

  rate = 1/a;

  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a)) z = R::rexp(rate);
  u = R::runif(0.0, 1.0);

  while( log(u) > (-0.5*z*z)){
    z = R::rexp(rate);
    while(z > (b-a)) z = R::rexp(rate);
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}


// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================
double rnorm_trunc (double mu, double sigma, double lower, double upper){
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;

  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;

  // First scenario
  if( (a == R_NegInf) || (b == R_PosInf)){
    if(a == R_NegInf){
      change = 1;
      a = -b;
      b = R_PosInf;
    }

    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  // Second scenario
  else if((a * b) <= 0.0){
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1)){
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  // Third scenario
  else{
    if(b < 0){
      tmp = b; b = -a; a = -tmp; change = 1;
    }

    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2) z = unif_rs(a,b);
    else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
    else z = exp_rs(a,b);
    if(change) z = -z;
  }
  double output;
  output = sigma*z + mu;
  return (output);
}


// Sample phi
// ======================================================================
// [[Rcpp::export]]
SEXP phiSampler(mat phi, mat theta, mat gamma, vec kappa, mat M,
                 mat Sigma, mat Sigmachol, vec nu){

  int p  = phi.n_cols;
  int q  = phi.n_rows;
  int np = gamma.n_rows;

  for(int j=0; j<p; j++ ){
    //Update theta_j
    for(int i=0; i<np; i++ ){
      double nu2   = 1/kappa(j) * nu(i);
      nu2          = pow(nu2, 0.5);
      double lower = gamma(i,j)*0 - (1-gamma(i,j))*1000000;
      double upper = gamma(i,j)*1000000 - (1-gamma(i,j))*0;
      theta(i,j)   = rnorm_trunc(0,nu2,lower,upper);
    }

    //Conditional on theta_j, update phi_j first
    vec phimean = Sigma*trans(M)*theta.col(j);
    vec normies = Rcpp::rnorm(q,0,1);
    phi.col(j) = phimean + sqrt(1/kappa(j))*Sigmachol*normies;
  }

  return Rcpp::List::create(Rcpp::Named("theta") = theta,
                            Rcpp::Named("phi")   = phi) ;
}



// Sample rho - the temporal correlation parameter
// ======================================================================
// [[Rcpp::export]]
mat rhoSampler(mat Y, mat X, mat beta, mat rho, vec sigma2){

  int T  = Y.n_rows;
  int P  = rho.n_cols;
  int n  = Y.n_cols;

  mat Id(P,P); Id.eye();
 
  for (int i=0; i<n; i++){
    vec Ystar = Y.col(i) - X*trans(beta.row(i));
    mat w((T-2),2); 
    w.col(0) = Ystar.rows(1, (T-2));
    w.col(1) = Ystar.rows(0, (T-3));

    mat Ginv    = inv(trans(w)*w + 3*Id);
    vec rhohat  = Ginv*trans(w)*Ystar.rows(2, (T-1));
    vec normies = Rcpp::rnorm(P,0,1);
    vec rho_    = rhohat + sqrt(sigma2(i))*trans(chol(Ginv))*normies;

    rho.row(i) = trans(rho_);
  }

  return(rho);
}


// Sample tau2 - spike/slab (slab) prior variance
// ======================================================================
// [[Rcpp::export]]
mat tau2Sampler(mat X, mat beta, mat gamma, vec tau2){

  int p  = X.n_cols;
  
  for (int j=0; j<p; j++){
    double t2shape = sum( gamma.col(j) )/2;
    vec t21        = gamma.col(j)%(beta.col(j)%beta.col(j));
    double t2scale = sum( t21 ) /2 + 0.001 ;
    double t2temp  = R::rgamma(t2shape+1.0, 1.0/t2scale);
    tau2(j)        = 1.0/t2temp;
  }

  return(tau2);
}


// Sample gamma and beta
// ======================================================================
// [[Rcpp::export]]
SEXP gammaSampler(mat Y, mat X, mat eta, mat beta, mat gamma, mat rho, 
                 vec tau2, vec sigma2){

  int n  = Y.n_cols;
  int p  = X.n_cols;
  int T  = Y.n_rows;
   
  for (int i=0; i<n; i++){
    for (int j=0; j<p; j++){
      //Update gamma
      vec Xj  = X.col(j); 
      vec XjL = Xj.rows(2,T-1) - Xj.rows(1,T-2)*rho(i,0) - Xj.rows(0,T-3)*rho(i,1);
      vec Yj  = Y.col(i) - X*trans(beta.row(i)) + Xj*beta(i,j);
      vec wj  = Yj.rows(2,T-1) - Yj.rows(1,T-2)*rho(i,0) - Yj.rows(0,T-3)*rho(i,1);

      double ee  = tau2(j) * sum(XjL%XjL) / sigma2(i) + 1;
      double BF1 = sqrt(ee);
      double e1  = sum(XjL%wj);
      double BF2 = exp( -0.5 * ( e1*e1/(sum(XjL%XjL) + sigma2(i)/tau2(j) )   ) / sigma2(i)  );
      double BF  = log(BF1) + log(BF2);

      double p1 = eta(i,j);
      double p0 = 1 - p1;

      gamma(i,j) = ::Rf_rbinom(1,p1/(p1+exp(BF)*p0) );

      //Update beta conditional on new gamma
      if (gamma(i,j) == 1 ){
        double betacov   = 1/( sum(XjL%XjL) + sigma2(i)/tau2(j)  ) ;            
        double betamean  = betacov * sum(XjL%wj);
        double betachol  = sqrt(sigma2(i)*betacov);
        //double normie    = ::Rf_rnorm(0,betachol);
        //beta(i,j)        = betamean + normie;
      }
      else{
        beta(i,j)        = 0;
      }

    }
  }
  
  return Rcpp::List::create(Rcpp::Named("gamma")  = gamma,
                            Rcpp::Named("beta")   = beta);
}



// Sample sigma2 - sample the error variance
// ======================================================================
// [[Rcpp::export]]
mat sigma2Sampler(mat Y, mat X, mat beta, mat rho){

  double n  = Y.n_cols;
  double T  = Y.n_rows;

  double gamA = T/2 - 1;

  vec sigma2(n);

  beta = trans(beta);

  for (int i=0; i<n; i++){
    vec Ystar = Y.col(i) - X*(beta.col(i));
    
    vec e1 = Ystar.rows(2,T-1) - Ystar.rows(1,T-2)*rho(i,0) - Ystar.rows(0,T-3)*rho(i,1);

    double temp = R::rgamma(gamA, 1.0/(0.5*(dot(e1,e1)) ));
    sigma2(i) = 1.0/temp;
  }

  return(sigma2);
}



// Sample kappa - spatial precision parameter
// ======================================================================
// [[Rcpp::export]]
mat kappaSampler(mat theta, vec nu){
 
  int p = theta.n_cols;
  int n = theta.n_rows;

  vec kappa = zeros<vec>(p);

  for (int j=0; j<p; j++){ 
    vec theta2  = theta.col(j) % theta.col(j);
    double temp = sum( theta2/nu);
    kappa(j)    = R::rgamma((n+1)/2, 1/(0.5*temp + 1/2000));
  }

  return(kappa);
}
