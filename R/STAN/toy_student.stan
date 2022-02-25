

// DATA
// 
data {
  int<lower=0> N; 
  vector[N] X;    
}

// PARAMETERS
// 
parameters {             
  real<lower=0> df;
  real mu;             
  real<lower=0> sigma;
}

// MODEL
//       
model {
  X ~ student_t(df, mu, sigma); 
  
  mu ~ normal(0, 1); 
  sigma ~ lognormal(0, 1);   
  df ~ lognormal(0, 1);   
}

// GENERATED QUANTITIES
//
generated quantities {
  vector[ N] logLikelihood; 
  
  for(i in 1:N) { 
    logLikelihood[i] = student_t_lpdf(X[i] | df, mu, sigma);
    }
}
