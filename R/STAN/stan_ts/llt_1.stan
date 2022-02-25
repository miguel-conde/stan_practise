// https://marubon-ds.blogspot.com/2018/01/local-linear-trend-model-for-time.html
// 
data {
  int N;
  vector[N] X;
}

parameters {
  vector[N] u;
  vector[N] v;
  real<lower=0> s_u;
  real<lower=0> s_v;
  real<lower=0> s_x;
}

model {
  v[2:N] ~ normal(v[1:N-1], s_v);
  u[2:N] ~ normal(u[1:N-1] + v[1:N-1], s_u);
  X ~ normal(u, s_x);
  
  s_u ~ normal(0, 1);
  s_v ~ normal(0, 1);
  s_x ~ normal(0, 1);
}

generated quantities {
  real logLikelihood[N];
  
  for(i in 1:N) {
    logLikelihood[i] = normal_lpdf(X[i] | u[i], s_x);
  }
}