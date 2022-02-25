
// FUNCTIONS
//
functions {
  real covariateMean(real aX, real aBeta) {
    return(aBeta * log(aX));
  }
}

// DATA
//
data {
  int<lower=0> N; 
  vector[N] y;    // Heights for N people 
  vector[N] x;    // weights
}

// PARAMETERS
//
parameters {
  real beta;            
  real<lower=0> sigma; 
}

// MODEL
//
model {
  
  for (i in 1:N) {
    y[i] ~ normal(covariateMean(x[i], beta), sigma); 
  }
  
  beta  ~ normal(0, 1); 
  sigma ~ gamma(1, 1);   
}

// GENERATED QUANTITIES
//
generated quantities {
  vector[N] lSimData;
  int aMax_indicator;
  int aMin_indicator;
  
  // Generate posterior predictive samples
  for (i in 1:N) {
    lSimData[i] = normal_rng(covariateMean(x[i], beta), sigma);
  }
  
  // Compare with real data
  aMax_indicator = max(lSimData) > max(y);
  aMin_indicator = min(lSimData) < max(y);
}


