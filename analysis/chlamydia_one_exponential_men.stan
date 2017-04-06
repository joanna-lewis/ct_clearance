// Fit a single exponential model to published data, to obtain estiamte of duration of untreated chlamydia infection in men

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
  real<lower=0> lambda; // clearance rate
  real<lower=0,upper=1> psi; // sensitivity of culture given initial positive culture
}

transformed parameters {
	
  real theta[Nobs]; // clearance probability
  
  // proportion of participants in each study time point expected to have recovered.
  
  // theta calculated the same way for all studies
  for (i in 1:studnum) { // should find a way of using input data to state study categories
    for (j in 1:studobs[i]) {
    
      theta[cumobs[i]+j] <-  1 - exp(-lambda*t[cumobs[i]+j]);
      
      if(seind[cumobs[i]+j])
        theta[cumobs[i]+j] <- 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }
       
  }

model {

  // priors
  lambda ~ exponential(0.001); // clearance rate
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
   
	
//# left truncation
//w1 <- (p1 / lambda.C[1]) / (p1 / lambda.C[1] + (1 - p1) / lambda.C[2])

//# summary statistics
//dur <- 1 / lambda.C[2]

//# Predicted values for Forest plot
//for (i in 1:studnum) {
// for (j in 1:studobs[i]) {	
//  stud.lambda.Cexpect[i,j] <- -log(1 - theta[i,j]) / t[i,j]
//  stud.dur.expect[i,j] <- 1 / stud.lambda.Cexpect[i,j]	
//  }
// }
//}

  
  }

