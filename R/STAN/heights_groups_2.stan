
data {
  int N;         // number of samples across all studies
  int K;         // number of studies
  real Y[N];     // heights for all N people
  int groups[N]; // id of each observation
}

parameters {
  real mu[K];                // mean height in population
  real<lower = 0> sigma[K]; // sd of height pop distribution
}

model {
  for (i in 1:K) {
    Y[i] ~ normal(mu[groups[i]], sigma[groups[i]]);
  }
  
  mu ~ normal(1.5, 0.1);
  sigma ~ gamma(1, 1);
}
