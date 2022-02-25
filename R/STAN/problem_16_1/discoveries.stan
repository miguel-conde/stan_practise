// Problem 16.1 Discoveries data revisited 

data {
  int<lower = 0> N; // Number of years
  int<lower = 0> X[N]; // Nummber of great inventions / discoveries each year
}

parameters {
  real<lower = 0> lambda; // Rate of occurence of discoveries
}

model {
  X ~ poisson(lambda);      // Likelihood
  
  lambda ~ lognormal(2, 1); // Prior
}

generated quantities {
  int<lower = 0> XSim[N];
  
  for (i in 1:N) {
    XSim[i] = poisson_rng(lambda);
  }
}
