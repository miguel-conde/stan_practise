data {
  
  int<lower=0> N;       // Number of observations
  int<lower=0> S;       // Seasonality
  vector[N] y;          // Response
}

transformed data {
    // Hidden States Process Transformation Matrix
  matrix[2, 2] B = [[1, 1], [0, 1]]; 
  // Observations Process Transformation Matrix
  matrix[1, 2] Z = [[1, 0]];
}


parameters {
  // Hidden States Process 
  vector<lower=0>[2] sigma_Q_llt; // Process errors covariance matrix
  vector<lower=0>[2] sigma_Q_seas;
  real w_llt[N];                     // Process errors
  real w_seas[N];
  
  // Observations Process
  real<lower=0> sigma_R;
}

transformed parameters {


  // 1 - Local Level Trend Compponent (LLT)
  // Hidden States level and slope
  vector[2] x_llt[N];
  // LLT forecasting
  real y_hat_llt[N];
  
  // 2 - Seasonal component
  real x_seas[N];
  
  // x_llt initial state
  x_llt[1] =[y[1], 0]'; 
  y_hat_llt[1] = y[1];
  
  for (i in 2:N) {
    x_llt[i] = B * x_llt[i-1] + w_llt[i];
    y_hat_llt[i] = (Z * x_llt[i])[1];
  }
  
  // // x_seas initial states
  x_seas[1] = 0;
  for (i in 2:(S-1)) {
    x_seas[i] = -sum(x_seas[1:(i-1)]) + w_seas[i];
  }

  for (i in S:N) {
    x_seas[i] = -sum(x_seas[(i-S+1):(i-1)])  + w_seas[i];
  }
  
}

model {
  
  for (i in 1:N) {
    
    w_llt[i] ~ normal(0, sigma_Q_llt);
    w_seas[i] ~ normal(0, sigma_Q_seas);
    
    y[i] ~ normal(y_hat_llt[i] + x_seas[i], sigma_R);
  }
  
  // sigma_Q_llt ~ gamma(1e3*sqrt(variance(y)), 10);
  // sigma_Q_seas ~ gamma(1e3*sqrt(variance(y)), 10);
  sigma_Q_llt ~ inv_gamma(0.1, 1);
  sigma_Q_seas ~ inv_gamma(0.1, 1);

  sigma_R ~ gamma(0.1, 1);
}

generated quantities {
  
  // vector[2] logLikelihood_x[N];
  vector[N] logLikelihood_y;
  
  for (i in 1:N) {
    // logLikelihood_x[i] = multi_normal_lpdf(x[i] | mu_x[i], quad_form_diag(Omega[i], tau[i]));
    logLikelihood_y[i] = normal_lpdf(y[i] | y_hat_llt[i] + x_seas[i], sigma_R);
  }
}
