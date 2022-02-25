// https://discourse.mc-stan.org/t/bayesian-structural-time-series-modeling/2256/2
// 
data {
  int <lower=0> T;
  vector[T] y;
}

parameters {
  vector[T] u_err; //Slope innovation
  vector[T] v_err; //Level innovation
  real <lower=0> s_obs;
  real <lower=0> s_slope;
  real <lower=0> s_level;
}

transformed parameters {
  vector[T] u; //Level
  vector[T] v; //Slope
  u[1] = u_err[1];
  v[1] = v_err[1];
  for (t in 2:T) {
    u[t] = u[t-1] + v[t-1] + s_level * u_err[t];
    v[t] = v[t-1] + s_slope * v_err[t];
  }
}

model {
  s_obs ~ normal(0,1);
  s_slope ~ normal(0,1);
  s_level ~ normal(0,1);
  
  u_err ~ normal(0,1);
  v_err ~ normal(0,1);
  y ~ normal (u, s_obs);
}

generated quantities {
  real logLikelihood[T];
  
  for(i in 1:T) {
    logLikelihood[i] = normal_lpdf(y[i] | u[i], s_obs);
  }
}