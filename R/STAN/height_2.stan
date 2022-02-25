



// DATA
// 
data {
  int<lower=0> N; 
  vector[N] y;    
}

// The final code block is known as transformed data, and as its name suggests 
// it can be used to make transformations to the data you pass to Stan. These 
// transformations are carried out once at the beginning of the program, so it 
// usually does not significantly affect the efficiency of execution 
// As a silly example, imagine that instead of fitting a model to the height 
// data itself, we want our model to explain the squared deviation from the 
// sample mean. One way to do this is to carry out the data transformation in 
// Stan:
//
// TRANSFORMED DATA
//             Block executed once after the data block is executed
//
transformed data {
  vector[N] lSqDeviation;
  
  for(i in 1:N) {
    lSqDeviation[i] = (y[i] - mean(y))^2;
  }
}

// PARAMETERS
// 
parameters {
  real mu;             
  real<lower=0> sigmaSq; // variance of height pop distribution 
}

// There are occasions when we want to generate samples for transformations of 
// those parameters defined in the parameters block and, possibly, even sample 
// from these transformed parameters. 
// In our original heights example (no covariates), imagine that instead of 
// setting priors on the standard deviation  parameters – sigma – you wish to do 
// so on the variance. However, you also want to generate samples for sigma. 
// One way to do this is to use the transformed parameters block:
//
// TRANSFORMED PARAMETERS
//                        Block executed each time the log probability is 
//                        evaluated
//
transformed parameters {
  real sigma;
  
  sigma = sqrt(sigmaSq);
}

// MODEL
//       
model {
  y ~ normal(mu, sigma); 
  
  mu ~ normal(1.5, 0.1); 
  sigmaSq ~ gamma(5, 1);   
}

// GENERATED QUANTITIES
//
generated quantities {
  vector[N] lSimData;
  int aMax_indicator;
  int aMin_indicator;
  
  // Generate posterior predictive samples
  for (i in 1:N) {
    lSimData[i] = normal_rng(mu, sigma);
  }
  
  // Compare with real data
  aMax_indicator = max(lSimData) > max(y);
  aMin_indicator = min(lSimData) < max(y);
}
