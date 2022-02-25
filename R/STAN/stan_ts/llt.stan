
data {
  
  int<lower=0> N;       // Number of observations
  vector[N] y;          // Response
}

transformed data {
  
}


parameters {
  vector[N] slope_err; // Slope innovation
  vector[N] level_err; // Level innovation
  
  real<lower=0> epsilon;
  real<lower=0> eta_lvl;
  real<lower=0> eta_slp;
}

transformed parameters {
  
  // real pred[N];
  
  real level[N];
  real slope[N];
  
  // pred[1] = y[1];
  
  level[1] = eta_lvl;
  slope[1] = eta_slp;
  
  for(i in 2:N) {
    
    level[i] = level[i-1] + slope[i-1] + eta_lvl*level_err[i]; 
    
    slope[i] = slope[i-1] + eta_slp*slope_err[i];
    
    // pred[i] = level[i]; 
  }
}


model {
  
  // y[2:N] ~ normal(pred[2:N], epsilon);
  y ~ normal(level, epsilon);
  
  slope_err ~ normal(0, 1);
  level_err ~ normal(0, 1);
  
  epsilon ~ normal(0, 1);
  eta_lvl ~ normal(0, 1);
  eta_slp ~ normal(0, 1);
}

generated quantities {
  real logLikelihood[N];
  
  for(i in 1:N) {
    logLikelihood[i] = normal_lpdf(y[i] | level[i], epsilon);
  }
}
