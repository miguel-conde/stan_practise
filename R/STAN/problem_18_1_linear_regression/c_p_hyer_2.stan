// State level model, time trend

data {
  int K;                                  // Number of state
  int N;                                  // Number of observations
  real murder[N];                         // Response
  int<lower = 0, upper = 1> penalty[N];   // Dummy var
  real car[N];                            // Regressor
  int state[N];                          // The state of each observation
  int t[N];
}


parameters {
  real alpha[K];
  real delta[K];
  real beta_car[K];
  real beta_penalty[K];
  real<lower = 0> sigma[K];
  
  real alpha_top;
  real<lower = 0> sigma_alpha_top;
  real delta_top;
  real<lower = 0> sigma_delta_top;
  real beta_car_top;
  real<lower = 0> sigma_car_top;
  real beta_penalty_top;
  real<lower = 0> sigma_penalty_top;
  
}


model {
  // Likelihood
  // for (i in 1:N) {
  //   murder[i] ~ normal(alpha + beta_car * car[i] + beta_penalty * penalty[i], sigma);
  // }
  
  for (i in 1:N) {
    murder[i] ~ normal(alpha[state[i]] + 
                          delta[state[i]] * t[i] + 
                          beta_car[state[i]] * car[i] + 
                          beta_penalty[state[i]] * penalty[i], 
                       sigma[state[i]]);
  }
  
  // Priors
  beta_car ~ normal(beta_car_top, sigma_car_top);
  beta_penalty ~ normal(beta_penalty_top, sigma_penalty_top);
  alpha ~ normal(alpha_top, sigma_alpha_top);
  delta ~ normal(delta_top, sigma_delta_top);
  sigma ~ normal(0, 1);
  
  // Hyper priors
   alpha_top ~ normal(0, 1);
   sigma_alpha_top ~ normal(0, 1);
   delta_top ~ normal(0, 1);
   sigma_delta_top ~ normal(0, 1);
   beta_car_top ~ normal(0, 1);
   sigma_car_top ~ normal(0, 1);
   beta_penalty_top ~ normal(0, 1);
   sigma_penalty_top ~ normal(0, 1);
}

generated quantities{
  real alpha_average;
  real delta_average;
  real beta_car_average;
  real beta_penalty_average;
  
  real logLikelihood[N];
  
  alpha_average = normal_rng(alpha_top, sigma_alpha_top);
  delta_average = normal_rng(delta_top, sigma_delta_top);
  beta_car_average = normal_rng(beta_car_top, sigma_car_top);
  beta_penalty_average = normal_rng(beta_penalty_top, sigma_penalty_top);
  
  for(i in 1:N) {
    logLikelihood[i] = normal_lpdf(murder[i] | alpha[state[i]] + 
                          delta[state[i]] * t[i] + 
                          beta_car[state[i]] * car[i] + 
                          beta_penalty[state[i]] * penalty[i], 
                       sigma[state[i]]);
  }
  
}
