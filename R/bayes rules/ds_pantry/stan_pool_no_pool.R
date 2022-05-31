library(tidyverse)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(broom.mixed)
library(rstan)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

# Data --------------------------------------------------------------------


# Load data
data(cherry_blossom_sample, package = "bayesrules")
running <- cherry_blossom_sample

# Remove NAs
running <- running %>% 
  select(runner, age, net) %>% 
  na.omit()

ggplot(running, aes(age, net)) + 
  geom_point() 

ggplot(running, aes(net)) + 
  geom_density()

ggplot(running, aes(age, net)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(running$runner)

ggplot(running) + 
  geom_density(aes(net, color = runner))


# POOLED MODEL ------------------------------------------------------------


stan_m1_code <- "
data {
  int<lower = 1> n;
  vector[n] age;
  vector[n] net;
}

transformed data {
  real avg_net;
  real sd_net;
  real avg_age;
  real sd_age;
  
  real s_alpha;
  real s_beta;
  
  avg_net = mean(net);
  sd_net = sd(net);
  avg_age = mean(age);
  sd_age = sd(age);
  
  s_alpha = sqrt(pow(sd_net, 2) + pow(sd_net / sd_age * avg_age, 2));
  s_beta  = sd_net / sd_age;
}

parameters {
  real alpha_centered;
  real beta_centered;
  
  real sigma_centered;

}

transformed parameters {
  real alpha;
  real beta;
  real sigma;
  
  vector[n] mu;
  
  alpha = avg_net + s_alpha * alpha_centered;
  beta  = 0 + s_beta * beta_centered;
  sigma = 0 + sd_net * sigma_centered;
  
  mu = alpha + beta * age;

}

model {  

  // PRIORS
  alpha_centered ~ normal(0, 1);
  beta_centered ~ normal(0, 1);
  sigma_centered ~ exponential(1);

  // MODEL
  net ~ normal(mu, sigma);
}

generated quantities {

  vector[n] mu_pps;
  vector[n] log_lik;

// PPS
  for (i in 1:n) {
    
    mu_pps[i] = alpha + beta * age[i];
    
    log_lik[i] = normal_lpdf(net[i] | mu_pps[i], sigma);
  }

}
"

stan_m1_code <- "
data {
  int<lower = 1> n;
  vector[n] age;
  vector[n] net;
}

parameters {
  real alpha;
  real beta;
  
  real sigma;

}

transformed parameters {
  vector[n] mu;
  
  mu = alpha + beta * age;

}

model {  

  // PRIORS
  alpha ~ normal(90, 300);
  beta ~ normal(1, 5);
  sigma ~ exponential(1);

  // MODEL
  net ~ normal(mu, sigma);
}

generated quantities {

  vector[n] mu_pps;
  vector[n] log_lik;

// PPS
  for (i in 1:n) {
    
    mu_pps[i] = alpha + beta * age[i];
    
    log_lik[i] = normal_lpdf(net[i] | mu_pps[i], sigma);
  }

}
"

stan_m1 <- stan_model(model_code = stan_m1_code)

fit_stan_m1 <- sampling(stan_m1, 
                        iter = 10000, chains = 1, 
                        # control = list(adapt_delta = 0.99),
                        data = compose_data(running),
                        verbose = TRUE)

the_pars <- c("alpha", "beta", "sigma")
print(fit_stan_m1, pars = the_pars)
plot(fit_stan_m1, pars = the_pars)

chains_stan_m1 <- extract(fit_stan_m1, the_pars) %>% bind_cols()
mcmc_trace(chains_stan_m1)
mcmc_dens_overlay(chains_stan_m1)
mcmc_acf(chains_stan_m1)
neff_ratio(fit_stan_m1)
rhat(fit_stan_m1)

loo(fit_stan_m1)
