data {
  int<lower=0> n; // Número de observaciones
  
  vector[n] height;

  // row_vector[m] L; // Lower bounds
  // row_vector[m] U; // Upper bounds
  
  vector[n] weight; 
  vector[n] age;
  vector[n] male;
  
  int intcpt;
}

transformed data {
  real avg_height;
  real sd_height;
  real avg_weight;
  real sd_weight;
  real avg_age;
  real sd_age;
  real avg_male;
  real sd_male;
  
  avg_height = mean(height);
  sd_height  = sd(height);
  avg_weight = mean(weight);
  sd_weight  = sd(weight);
  avg_age    = mean(age);
  sd_age     = sd(age);
  avg_male   = mean(male);
  sd_male    = sd(male);
  
}


parameters {
  real scaled_beta_weight;
  real scaled_beta_age;
  real scaled_beta_male;
  real scaled_beta_intcpt;
  real<lower=0> scaled_sigma;
}

transformed parameters {
  real beta_weight;
  real beta_age;
  real beta_male;
  real beta_intcpt;
  real<lower=0> sigma;
  
  vector[n] mu;
  
  beta_weight = sd_height / sd_weight * scaled_beta_weight;
  beta_age    = sd_height / sd_age    * scaled_beta_age;
  beta_male   = sd_height / sd_male   * scaled_beta_male;
  if (intcpt == 1) {
    beta_intcpt = avg_height + sqrt(pow(sd_height, 2) +
                                    pow(sd_height/sd_weight*avg_weight, 2) +
                                    pow(sd_height/sd_age*avg_age, 2) +
                                    pow(sd_height/sd_male*avg_male, 2)) * scaled_beta_intcpt;
  } else {
    beta_intcpt = 0;
  }
  sigma = sd_height * scaled_sigma;
  
  mu = beta_intcpt + beta_weight * weight + beta_age * age + beta_male * male;
  
}


model {
  
  // PRIORS
  scaled_beta_weight ~ normal(0, 1);
  scaled_beta_age    ~ normal(0, 1);
  scaled_beta_male   ~ normal(0, 1);
  scaled_beta_intcpt ~ normal(0, 1);
  
  scaled_sigma       ~ exponential(1);
  
  
  // LIKELIHOOD
  height ~ normal(mu, sigma);
}

