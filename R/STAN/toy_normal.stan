

// DATA
// 
data {
  int<lower=0> N; 
  vector[N] X;    
}

// PARAMETERS
// 
parameters {
  real mu;             
  real<lower=0> sigma;
}

// MODEL
//       
model {
  X ~ normal(mu, sigma); 
  
  mu ~ normal(0, 1); 
  sigma ~ lognormal(0, 1);   
}

// GENERATED QUANTITIES
//
generated quantities {
  vector[ N] logLikelihood; 
  
  for(i in 1: N) { 
    logLikelihood[i] = normal_lpdf(X[i] | mu, sigma); 
    }
}
