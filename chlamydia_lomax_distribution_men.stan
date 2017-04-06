// Based on Supplementary Web Appendix 1 to:
// Malcolm J. Price, A. E. Ades, Daniela De Angelis, Nicky J. Welton, John Macleod, Kate Soldan, Katy Turner, Ian Simms and Paddy J. Horner
// Mixture-of-exponentials models to explain heterogeneity in studies of the duration of Chlamydia trachomatis infection.
// Statistics in Medicine 32:1547–1560 (2013)

// duration of infection in men - no left-truncated, repeat-observation studies

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
  real<lower=1> k; // shape for gamma distribution of rates
  real<lower=0> alpha; // scale for gamma distribution of rates
  real<lower=0,upper=1> psi; // sensitivity of culture given initial positive culture
}

transformed parameters {
	
  real theta[Nobs]; // weighted average of clearance probabilities for two classes
    
  // proportion of participants in each study time point expected to have recovered.
  // clinic studies
  for (i in 1:studnum_bytype[1]) { 
    for (j in 1:studobs[i]) {

			theta[cumobs[i]+j] = pareto_type_2_cdf(t[cumobs[i]+j], 0, 1/k, alpha);

            if(seind[cumobs[i]+j])
              theta[cumobs[i]+j] = 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }
   
  // left-truncated; single observation
  for (i in (studnum_bytype[1]+1):(studnum_bytype[1]+studnum_bytype[2])) { 
    for (j in 1:studobs[i]) {

			theta[cumobs[i]+j] = pareto_type_2_cdf(t[cumobs[i]+j], 0, 1/(k-1), alpha);

            if(seind[cumobs[i]+j])
              theta[cumobs[i]+j] = 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }  
    
  }

model {

  // priors
  k ~ exponential(0.001); // shape for gamma distribution of rates
  alpha ~ exponential(0.001); // scale for gamma distribution of rates
  psi ~ beta(78,8); // sensitivity of culture given initial positive culture
	
  // likelihood
  for (i in 1:studnum) {
    for (j in 1:studobs[i]) {
//      print(i);
//      print("--->", j);
      r[cumobs[i]+j] ~ binomial(n[cumobs[i]+j], theta[cumobs[i]+j]);
      }
    }

  }

generated quantities {

  real dev[Nobs];
  real sumdev;
  int<lower=0> r_sim[Nobs]; // simulated number who cleared infection at each time point 

  // deviance
  for (i in 1:studnum) {
    for (j in 1:studobs[i]) {
    
      if(r[cumobs[i]+j] == 0)
      	
      	dev[cumobs[i]+j] = 2 * ( (n[cumobs[i]+j] - r[cumobs[i]+j]) * log((n[cumobs[i]+j] - r[cumobs[i]+j]) / (n[cumobs[i]+j] - (n[cumobs[i]+j] * theta[cumobs[i]+j]))));
      	
      else if(r[cumobs[i]+j] == n[cumobs[i]+j])
      
      dev[cumobs[i]+j] = 2 * (r[cumobs[i]+j] * log(r[cumobs[i]+j] / (theta[cumobs[i]+j] * n[cumobs[i]+j])));
      	
      else	
   
        dev[cumobs[i]+j] = 2 * (r[cumobs[i]+j] * log(r[cumobs[i]+j] / (theta[cumobs[i]+j] * n[cumobs[i]+j])) + (n[cumobs[i]+j] - r[cumobs[i]+j]) * log((n[cumobs[i]+j] - r[cumobs[i]+j]) / (n[cumobs[i]+j] - (n[cumobs[i]+j] * theta[cumobs[i]+j]))));	
      
      }
       
     }
     
  sumdev = sum(dev);
  
  // simulated number clearing infection (for Forest plot)
  for (i in 1:studnum) {
    for (j in 1:studobs[i]) {
      r_sim[cumobs[i]+j] = binomial_rng(n[cumobs[i]+j], theta[cumobs[i]+j]);
      }
    }
  
  }

