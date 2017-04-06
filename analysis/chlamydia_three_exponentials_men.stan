// Based on Supplementary Web Appendix 1 to:
// Malcolm J. Price, A. E. Ades, Daniela De Angelis, Nicky J. Welton, John Macleod, Kate Soldan, Katy Turner, Ian Simms and Paddy J. Horner
// Mixture-of-exponentials models to explain heterogeneity in studies of the duration of Chlamydia trachomatis infection.
// Statistics in Medicine 32:1547–1560 (2013)

data {
  int<lower=0> studnum; // number of studies 
  int<lower=0> studnum_bytype[3]; // number of studies 
  int<lower=0> studobs[studnum]; // number of observations (time periods) in each study
  int<lower=0> cumobs[studnum]; // cumulative number of observations at the start of each study
  int<lower=0> Nobs; // total number of observations (=sum(studobs))
  int<lower=0> r[Nobs]; // number who cleared infection at each time point 
  int<lower=0> n[Nobs]; // number tested at each time point 
  real<lower=0> t[Nobs]; // estimated mean follow-up 
  int<lower=0> seind[Nobs]; // did the study use culture as opposed to NAATs?
  real<lower=0> T[Nobs]; // already followed for...
}


parameters {
  real<lower=0,upper=1> p1; // probability of clearing fast
  real<lower=0,upper=1> p2; // probability of clearing at medium rate
  real<lower=0> lambda_mid; // medium clearance rates
  real<lower=0> lambda_slow; // slow clearance rates
  real<lower=0,upper=1> psi; // sensitivity of culture given initial positive culture
}

transformed parameters {
	
  real<lower=0,upper=1> w1; //  expected proportion of fast clearers in a prevalent population
  real<lower=0,upper=1> w2; //  expected proportion of medium clearers in a prevalent population
  real theta[Nobs]; // weighted average of clearance probabilities for two classes
  real<lower=0> lambda[3];
  vector[3] p; // vector to contain prababilities of clearance rates
  real<lower=0,upper=1> pkm; // temporary variable, containing probability of this k, m
  
  real narray[3]; // temporary array, containing current combination of clearance class allocations
  
  // rates
  lambda[1] <- 49;
  lambda[2] <- lambda_mid;
  lambda[3] <- lambda_slow;
  
  // probabilities of each rate
  p[1] <- p1;
  p[2] <- p2;
  p[3] <- 1 - p1 - p2;
  
  // proportion of participants in each study time point expected to have recovered.
  
  // clinic studies
  
  for (i in 1:studnum_bytype[1]) { // should find a way of using input data to state study categories
    for (j in 1:studobs[i]) {
      theta[cumobs[i]+j] <-  0;
      for(k in 0:n[cumobs[i]+j]){ // k indexes fast clearers
        narray[1] <- k;
        for(m in 0:(n[cumobs[i]+j] - k)){ //m indexed medium clearers, from those who don't clear fast
          narray[2] <- m;
          narray[3] <- n[cumobs[i]+j]-k-m;
      
//          if(k == 0)
          
//            pkm <- (falling_factorial(n[cumobs[i]+j],n[cumobs[i]+j])/(falling_factorial(m,m)*falling_factorial(n[cumobs[i]+j]-m,n[cumobs[i]+j]-m)) ) * p[2]^m * p[3]^(n[cumobs[i]+j] - m);
            
//          else if(m == 0)
      
//            pkm <- (falling_factorial(n[cumobs[i]+j],n[cumobs[i]+j])/(falling_factorial(k,k)*falling_factorial(n[cumobs[i]+j]-k,n[cumobs[i]+j]-k)) )* p[1]^k * p[3]^(n[cumobs[i]+j] - k);
          
//          else if((k + m) == n[cumobs[i]+j])
          
//            pkm <-  (falling_factorial(n[cumobs[i]+j],n[cumobs[i]+j])/(falling_factorial(k,k)*falling_factorial(m,m)) )* p[1]^k * p[2]^m;
          
//          else{
            pkm <-  (falling_factorial(n[cumobs[i]+j],n[cumobs[i]+j])/(falling_factorial(k,k)*falling_factorial(m,m)*falling_factorial(n[cumobs[i]+j]-k-m,n[cumobs[i]+j]-k-m)) )* p[1]^k * p[2]^m * p[3]^(n[cumobs[i]+j] - k - m);
//            }
      
          theta[cumobs[i]+j] <- theta[cumobs[i]+j] + pkm * ( k/(1.0*n[cumobs[i]+j]) * (1 - exp(-lambda[1]*t[cumobs[i]+j])) + m/(1.0*n[cumobs[i]+j]) * (1 - exp(-lambda[2]*t[cumobs[i]+j])) + (n[cumobs[i]+j] - k - m)/(1.0*n[cumobs[i]+j]) * (1 - exp(-lambda[3]*t[cumobs[i]+j])) );
        
          }
        }
      if(seind[cumobs[i]+j])
        theta[cumobs[i]+j] <- 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }
   
  // left-truncated; single observation
  
  w1 <- (p1 / lambda[1]) / (p1 / lambda[1] + p2 / lambda[2] + (1-p1-p2) / (lambda[3]));
  w2 <- (p2 / lambda[1]) / (p1 / lambda[1] + p2 / lambda[2] + (1-p1-p2) / (lambda[3]));
  
  
  for (i in (studnum_bytype[1]+1):(studnum_bytype[1]+studnum_bytype[2])) { // should find a way of using input data to state study categories
    for (j in 1:studobs[i]) {
      theta[cumobs[i]+j] <-  0;
      for(k in 0:n[cumobs[i]+j]){ // k indexes fast clearers
        for(m in 0:(n[cumobs[i]+j]-k)){ //m indexed medium clearers, from those who don't clear fast
      
          pkm <- (falling_factorial(n[cumobs[i]+j],n[cumobs[i]+j])/(falling_factorial(k,k)*falling_factorial(m,m)*falling_factorial(n[cumobs[i]+j]-k-m,n[cumobs[i]+j]-k-m)) )* p[1]^k * p[2]^m * p[3]^(n[cumobs[i]+j] - k - m);
      
        theta[cumobs[i]+j] <- theta[cumobs[i]+j] + pkm * ( k/(1.0*n[cumobs[i]+j]) * (1 - exp(-lambda[1]*t[cumobs[i]+j])) + m/(1.0*n[cumobs[i]+j]) * (1 - exp(-lambda[2]*t[cumobs[i]+j])) + (n[cumobs[i]+j] - k - m)/(1.0*n[cumobs[i]+j]) * (1 - exp(-lambda[3]*t[cumobs[i]+j])) );
        
          }
        }
      if(seind[cumobs[i]+j])
        theta[cumobs[i]+j] <- 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }
      
    
  }

model {

  vector[3] alpha; // for diriclet prior

  alpha[1] <- 1;
  alpha[2] <- 1;
  alpha[3] <- 1;

  // priors
  p ~ dirichlet(alpha);
//  increment_log_prob(uniform_log(log(lambda_mid), -1, 10)); // medium clearance rate
  increment_log_prob(uniform_log(log(lambda_slow), -2.62, log(lambda_mid))); // slow clearance rate
  psi ~ beta(78,8); // sensitivity of culture given initial positive culture	

  // likelihood
  for (i in 1:studnum) {
    for (j in 1:studobs[i]) {
//      print(i);
//      print(j);
      r[cumobs[i]+j] ~ binomial(n[cumobs[i]+j], theta[cumobs[i]+j]);
      }
    }

  }

generated quantities {


  real dev[Nobs];
  real sumdev;

  // deviance
  for (i in 1:studnum) {
    for (j in 1:studobs[i]) {
    
      if(r[cumobs[i]+j] == 0)
      	
      	dev[cumobs[i]+j] <- 2 * ( (n[cumobs[i]+j] - r[cumobs[i]+j]) * log((n[cumobs[i]+j] - r[cumobs[i]+j]) / (n[cumobs[i]+j] - (n[cumobs[i]+j] * theta[cumobs[i]+j]))));
      	
      else if(r[cumobs[i]+j] == n[cumobs[i]+j])
      
      dev[cumobs[i]+j] <- 2 * (r[cumobs[i]+j] * log(r[cumobs[i]+j] / (theta[cumobs[i]+j] * n[cumobs[i]+j])));
      	
      else	
   
        dev[cumobs[i]+j] <- 2 * (r[cumobs[i]+j] * log(r[cumobs[i]+j] / (theta[cumobs[i]+j] * n[cumobs[i]+j])) + (n[cumobs[i]+j] - r[cumobs[i]+j]) * log((n[cumobs[i]+j] - r[cumobs[i]+j]) / (n[cumobs[i]+j] - (n[cumobs[i]+j] * theta[cumobs[i]+j]))));	
      
      }
       
     }
     
  sumdev <- sum(dev);
     
}

