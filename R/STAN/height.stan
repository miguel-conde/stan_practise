//
// Suppose that we want to infer the mean height, µ, in a population of interest.
//
// Suppose that the normal sampling model is appropriate given the multitude of 
// different factors that affect the height of a person. σ is the standard 
// deviation of the sampling distribution.
//
// To complete our specification of a Bayesian model we need priors, which we 
// choose as: 

// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// DATA
//      Block executed once at the beginning
//      Data types: real, int, vector, matrix, array
//                  (can be BOUNDED)
//
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // A discrete variable that takes integer values with a 
                  // minimum value of 0
  vector[N] y;    // Heights for N people - a vector of continuous variables of 
                  //length N
}

// PARAMETERS: 
//             Block executed each time the log probability is evaluated
//             More exotic data types are available here: simplex, corr_matrix,
//             ordered (vector), ...
//
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;             // Mean height in population - a continuous variable with  
                       // a minimum value of 0
  real<lower=0> sigma; // sd of height distribution - a continuous variable
                       //
}

// MODEL
//       Block executed each time the log probability is evaluated
//       Used to specify the likelihood and priors for a model
//
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // Y ~ normal( mu, sigma), which means ‘increment the overall log probability 
  // by an amount given by the log likelihood of a data point Y for a normal 
  // distribution with a mean of mu and standard deviation of sigma’
  y ~ normal(mu, sigma); // vectorised sampling statement for likelihood
  
  mu ~ normal(1.5, 0.1); // vectorised sampling statement for prior for mu
  sigma ~ gamma(1, 1);   // vectorised sampling statement for prior for sigma
}

// GENERATED QUANTITIES
//                      Block executed once per sample
//
// One of the main uses of this section is to do posterior predictive checks of 
// a model’s fit. Once you know how this section works, you will find that it is 
// much easier to do posterior predictive checks here rather than afterwards in 
// R, Python, and so on. 
// This section of a Stan program is executed once per sample, meaning that it 
// does not typically pose a threat to efficiency (although it can significantly 
// affect the memory used by Stan.
//
// Another use of this code block is to generate the requisite data for 
// measuring a model’s predictive performance.
//
// Finally, another use is that the generated quantities block can be used to 
// generate samples from parameters  that interest us at a given level of a 
// hierarchical model
//                      
// Let’s use our heights example from the preceding sections to illustrate how 
// we can use this code block to do a posterior predictive check. Suppose that 
// we want to test whether our model can generate the extremes that we see in 
// our data. In particular, we choose to count the fraction of posterior 
// predictive samples – each of the same size as our original data – where the 
// maximum or minimum of the simulated data is more extreme than the actual data.
generated quantities {
  vector[N] lSimData;
  int aMax_indicator;
  int aMin_indicator;
  
  // Generate posterior predictive samples
  for (i in 1:N) {
    // Note that to generate random samples from a given distribution, we use 
    // the ‘_rng’ suffix.
    // Y = normal_rng( mu, sigma) generates a single (pseudo-) independent 
    // sample from a normal distribution with a mean of mu and a standard 
    // deviation of sigma.
    lSimData[i] = normal_rng(mu, sigma);
  }
  
  // Compare with real data
  aMax_indicator = max(lSimData) > max(y);
  aMin_indicator = min(lSimData) < max(y);
}

