// Homogeneous model

data {
  int N;                                  // Number of observations
  real murder[N];                         // Response
  int<lower = 0, upper = 1> penalty[N];   // Dummy var
  real car[N];                            // Regressor
}


parameters {
  real alpha;
  real beta_car;
  real beta_penalty;
  real<lower = 0> sigma;
}


model {
  // Likelihood
  for (i in 1:N) {
    murder[i] ~ normal(alpha + beta_car * car[i] + beta_penalty * penalty[i], sigma);
  }
  
  // Priors
  beta_car ~ normal(0, 1);
  beta_penalty ~ normal(0, 1);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
}

generated quantities {
  real logLikelihood[N];
  
  for(i in 1:N) {
    logLikelihood[i] = normal_lpdf(murder[i] | alpha + beta_penalty * penalty[i] +
beta_car * car[i], sigma);
  }
}
