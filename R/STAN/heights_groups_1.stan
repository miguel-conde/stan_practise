
data {
  int N;        // number of samples across all studies
  int K;        // number of studies
  real Y[N];    // heights for all N people
  int S[K];     // sample sizes of each study
  int index[K]; // start position of each study
}

parameters {
  real mu[K];                // mean height in population
  real<lower = 0> sigma[K]; // sd of height pop distribution
}

model {
  for (i in 1:K) {
    Y[index[i]:index[i]+S[i]-1] ~ normal(mu[i], sigma[i]);
  }
  
  mu ~ normal(1.5, 0.1);
  sigma ~ gamma(1, 1);
}

