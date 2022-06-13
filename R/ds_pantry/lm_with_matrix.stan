data {
  int<lower=0> n; // Número de observaciones
  vector[n] y;
  
  int<lower=1> m; // Número de predictores
  // row_vector[m] L; // Lower bounds
  // row_vector[m] U; // Upper bounds
  
  row_vector[m] X[n]; // m "columnas" (predictores) x n "filas" (observaciones)
}


parameters {
  row_vector[m] betas; 
  real<lower=0> sigma;
}

transformed parameters {
  vector[n] mu;
  
  for ( i in 1:n) {
    mu[i] = dot_product(X[i], betas);
  }
  
}


model {
  
  // PRIORS
  sigma ~ exponential(10);
  
  for (i in 1:m) {
    betas[i] ~ normal(0, 10);
  }
  
  
  // LIKELIHOOD
  y ~ normal(mu, sigma);
}

