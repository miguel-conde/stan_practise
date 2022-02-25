//
// Sometimes it is helpful to use Stan as a tool to generate independent samples
// from a distribution of interest. This may be useful for posterior predictive 
// checking or, alternatively, because we want to know how a given Stan 
// distribution behaves. Suppose that we want to know what independent samples 
// from the Stan neg_binomial_2 distribution look like.

data {
  real mu;
  real kappa;
}

model {
  
}

generated quantities {
  int y;
  
  y = neg_binomial_2_rng(mu, kappa);
}
