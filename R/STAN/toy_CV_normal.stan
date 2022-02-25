
// DATA
// 
data {
  int<lower=0> NTest; 
  vector[NTest] XTest;   
  int<lower=0> NTrain; 
  vector[NTrain] XTrain;
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
  XTrain ~ normal(mu, sigma); 
  
  mu ~ normal(0, 1); 
  sigma ~ lognormal(0, 1);   
}

// GENERATED QUANTITIES
//
generated quantities {
  vector[NTest] logLikelihood; 
  
  for(i in 1:NTest) { 
    logLikelihood[i] = normal_lpdf(XTest[i] | mu, sigma); 
    }
}
