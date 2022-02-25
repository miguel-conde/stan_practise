
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
  real<lower=0> df;
  real mu;             
  real<lower=0> sigma;
}

// MODEL
//       
model {
  XTrain ~ student_t(df, mu, sigma); 
  
  mu ~ normal(0, 1); 
  sigma ~ lognormal(0, 1);   
  df ~ lognormal(0, 1);     
}

// GENERATED QUANTITIES
//
generated quantities {
  vector[NTest] logLikelihood; 
  
  for(i in 1:NTest) { 
    logLikelihood[i] = student_t_lpdf(XTest[i] | df, mu, sigma);
    }
}

