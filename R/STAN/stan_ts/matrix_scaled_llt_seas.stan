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
  real<lower = 0> std_dev_y;
  
  vector[N] y_scaled;
  
  std_dev_y = sd(y);
  y_scaled = y / std_dev_y;
  
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
  vector[2] x_scaled_llt[N];
  vector[2] x_llt[N];
  // LLT forecasting
  real y_hat_scaled_llt[N];
  real y_hat_llt[N];
  
  // 2 - Seasonal component
  real x_seas_scaled[N];
  real x_seas[N];
  
  real mu_scaled[N];
  real mu[N];
  
  // x_llt initial state
  x_scaled_llt[1] = [y_scaled[1], 0]'; 
  x_llt[1] = std_dev_y * x_scaled_llt[1];
  y_hat_scaled_llt[1] = y_scaled[1];
  y_hat_llt[1] = std_dev_y * y_hat_scaled_llt[1];
  
  // x_seas initial states
  x_seas_scaled[1] = 0;
  x_seas[1] = 0;
  
  mu_scaled[1] = y_hat_scaled_llt[1] + x_seas_scaled[1];
  mu[1] = std_dev_y * mu_scaled[1];
  
  for (i in 2:(S-1)) {
    
    x_scaled_llt[i] = B * x_scaled_llt[i-1] + w_llt[i];
    x_llt[i] = std_dev_y * x_scaled_llt[i];
    y_hat_scaled_llt[i] = (Z * x_scaled_llt[i])[1];
    y_hat_llt[i] = std_dev_y * y_hat_scaled_llt[i];
    
    x_seas_scaled[i] = -sum(x_seas_scaled[1:(i-1)]) + w_seas[i];
    x_seas[i] = std_dev_y * x_seas_scaled[i];
    
    mu_scaled[i] = y_hat_scaled_llt[i] + x_seas_scaled[i];
    mu[i] = std_dev_y * mu_scaled[i];
  }

  for (i in S:N) {
    
    x_scaled_llt[i] = B * x_scaled_llt[i-1] + w_llt[i];
    x_llt[i] = std_dev_y * x_scaled_llt[i];
    y_hat_scaled_llt[i] = (Z * x_scaled_llt[i])[1];
    y_hat_llt[i] = std_dev_y * y_hat_scaled_llt[i];
    
    x_seas_scaled[i] = -sum(x_seas_scaled[(i-S+1):(i-1)])  + w_seas[i];
    x_seas[i] = std_dev_y * x_seas_scaled[i];
    
    mu_scaled[i] = y_hat_scaled_llt[i] + x_seas_scaled[i];
    mu[i] = std_dev_y * mu_scaled[i];
  }
}

model {
  
  for (i in 1:N) {
    
    w_llt[i] ~ normal(0, sigma_Q_llt);
    w_seas[i] ~ normal(0, sigma_Q_seas);
    
    y[i] ~ normal(mu[i], sigma_R);
  }
  
  sigma_Q_llt ~ inv_gamma(1, 1);
  sigma_Q_seas ~ inv_gamma(1, 1);

  sigma_R ~ gamma(0.1, 1);
}

generated quantities {
  
  vector[N] logLikelihood_y;
  
  for (i in 1:N) {
    logLikelihood_y[i] = normal_lpdf(y[i] | mu[i], sigma_R);
  }
}

