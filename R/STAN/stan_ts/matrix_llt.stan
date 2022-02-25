data {
  
  int<lower=0> N;       // Number of observations
  vector[N] y;          // Response
}

transformed data {
  
}


parameters {
  vector<lower=0>[2] sigma_Q[N];
  real<lower=0> sigma_R[N];
  real w[N];
}

transformed parameters {
  matrix[2, 2] B = [[1, 1], [0, 1]];
  matrix[1, 2] Z = [[1, 0]];

  vector[2] x[N];
  
  x[1] =[y[1], 0]'; 
  
  for (i in 2:N) {
    x[i] = B * x[i-1] + w[i];
  }
}

model {
  
  for (i in 1:N) {
    
    sigma_Q[i] ~ gamma(1e3*sqrt(variance(y)), 10);
    // sigma_Q[i][1] ~ gamma(0.1, 1);
    // sigma_Q[i][2] ~ gamma(0.1, 1);
    
    w[i] ~ normal(0, sigma_Q[i]);
    
    sigma_R[i] ~ gamma(0.1, 1);
    
    y[i] ~ normal((Z * x[i])[1], sigma_R[i]);
  }
}

generated quantities {
  
  // vector[2] logLikelihood_x[N];
  vector[N] logLikelihood_y;
  
  for (i in 1:N) {
    // logLikelihood_x[i] = multi_normal_lpdf(x[i] | mu_x[i], quad_form_diag(Omega[i], tau[i]));
    logLikelihood_y[i] = normal_lpdf(y[i] | (Z * x[i])[1], sigma_R[i]);
  }
}
