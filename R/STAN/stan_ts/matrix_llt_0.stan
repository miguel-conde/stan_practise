data {
  
  int<lower=0> N;       // Number of observations
  vector[N] y;          // Response
}

transformed data {
  
}


parameters {
  // cov_matrix[2] Q[N];
  real<lower=0> R[N];
  vector<lower=0>[2] tau[N];
  corr_matrix[2] Omega[N];
  
  vector[2] x[N]; // Hidden states
}

transformed parameters {
  matrix[2, 2] B = [[1, 1], [0, 1]];
  matrix[1, 2] Z = [[1, 0]];

  vector[2] mu_x[N];
  real mu_y[N];
  
  mu_x[1] =[y[1], 0]'; 
  mu_y[1] = y[1];
  
  for (i in 2:N) {
    mu_x[i] = B * x[i-1];
    mu_y[i] = (Z * x[i])[1];
  }
}

model {
  
  for (i in 1:N) {
    
    // Q[i] ~ wishart(200, [[100, 0], [0, 100]]);
    tau[i] ~ cauchy(0, 2.5);
    Omega[i] ~ lkj_corr(2);
    
    R[i] ~ gamma(0.1, 1);
    
    // x[i] ~ multi_normal(mu_x[i]', sqrt(Q[i]));
    x[i] ~ multi_normal(mu_x[i]', quad_form_diag(Omega[i], tau[i]));
    
    y[i] ~ normal(mu_y[i], R[i]);
  }
}

generated quantities {
  
  // vector[2] logLikelihood_x[N];
  vector[N] logLikelihood_y;
  
  for (i in 1:N) {
    // logLikelihood_x[i] = multi_normal_lpdf(x[i] | mu_x[i], quad_form_diag(Omega[i], tau[i]));
    logLikelihood_y[i] = normal_lpdf(y[i] | mu_y[i], R[i]);
  }
}
