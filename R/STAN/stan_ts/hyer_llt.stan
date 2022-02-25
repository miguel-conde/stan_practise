data {
  
  int<lower=0> n;       // Number of observations
  int<lower=0> n_country;
  vector[n] y;          // Response
  int country[n];
  
}

transformed data {
  
  int<lower=0> n_obs_ctry[n_country];
  
  for(k in 1:n_country) {
    n_obs_ctry[k] = 0;
  }
  
  for(i in 1:n) {
    n_obs_ctry[country[i]] += 1;
  }
  
  print(n_obs_ctry);
}


parameters {
  real slope_err[n_country, n]; // Slope innovation
  real level_err[n_country, n]; // Level innovation
  
  real<lower=0> epsilon[n_country];
  real<lower=0> eta_lvl[n_country];
  real<lower=0> eta_slp[n_country];
  
  real slope_err_top; 
  real<lower = 0> sigma_slope_err_top;
  real level_err_top; 
  real<lower = 0> sigma_level_err_top;
  
  real epsilon_top;
  real<lower = 0> sigma_epsilon_top;
  real eta_lvl_top;
  real<lower = 0> sigma_eta_lvl_top;
  real eta_slp_top;
  real<lower = 0> sigma_eta_slp_top;
}

transformed parameters {
  
  vector[n]  level[n_country];
  vector[n]  slope[n_country];
  
  for(ctry in 1:n_country) {
    level[ctry][1] = eta_lvl[ctry];
    slope[ctry][1] = eta_slp[ctry];
    
    for(i in 2:n_obs_ctry[ctry]) {
      
      level[ctry][i] = level[ctry][i-1] + slope[ctry][i-1] + eta_lvl[ctry]*level_err[ctry][i]; 
      
      slope[ctry][i] = slope[ctry][i-1] + eta_slp[ctry]*slope_err[ctry][i];
    }
  }
}


model {
  
  for(i in 1:n) {
    
    y[i] ~ normal(level[country[i]], epsilon[country[i]]);
    
    slope_err[country[i]] ~ normal(slope_err_top, sigma_slope_err_top);
    level_err[country[i]] ~ normal(level_err_top, sigma_level_err_top);
    
    epsilon[country[i]] ~ normal(epsilon_top, sigma_epsilon_top);
    eta_lvl[country[i]] ~ normal(eta_lvl_top, sigma_eta_lvl_top);
    eta_slp[country[i]] ~ normal(eta_slp_top, sigma_eta_slp_top);
  }
  
  slope_err_top ~ normal(0, 1); 
  sigma_slope_err_top ~ gamma(1, 1);
  level_err_top ~ normal(0, 1); 
  sigma_level_err_top ~ gamma(1, 1);
  
  epsilon_top ~ normal(0, 1);
  sigma_epsilon_top ~ gamma(1, 1);
  eta_lvl_top ~ normal(0, 1);
  sigma_eta_lvl_top ~ gamma(1, 1);
  eta_slp_top ~ normal(0, 1);
  sigma_eta_slp_top ~ gamma(1, 1);
}

generated quantities {

  real<lower=0> epsilon_average;
  real<lower=0> eta_lvl_average;
  real<lower=0> eta_slp_average;

  real logLikelihood[n];

  epsilon_average = normal_rng(epsilon_top, sigma_epsilon_top);
  eta_lvl_average = normal_rng(eta_lvl_top, sigma_eta_lvl_top);
  eta_slp_average = normal_rng(eta_slp_top, sigma_eta_slp_top);

  for(i in 1:n) {
    logLikelihood[i] = normal_lpdf(y[i] | level[i], epsilon);
  }
}
