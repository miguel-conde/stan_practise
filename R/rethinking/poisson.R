library(tidyverse)
library(rethinking)
library(brms)
library(tidybayes)
library(bayesplot)
library(rstan)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

data(Kline)
d <- Kline %>% as_data_frame()
d

d <- d %>% 
  mutate(P = scale(log(population))[,1],
         contact_id = ifelse(contact == "high", 2, 1))
d


poisson_m1_code <- "
data {

  int<lower = 1> n;
  int<lower = 1> N_contacts;
  int cid[n];                // contact id
  real P[n];                 // Scaled log-Population
  
  int<lower = 0> total_tools[n];
}

parameters {

  vector[N_contacts] alpha;
  vector[N_contacts] beta;
}

transformed parameters {

  real log_lambda[n];
  real<lower = 0> lambda[n];
  
  for (i in 1:n) {
    log_lambda[i] = alpha[cid[i]] + beta[cid[i]] * P[i];
    lambda[i] = exp(log_lambda[i]);
  }
}

model {

  alpha ~ normal(3, 0.5);
  beta ~ normal(0, 0.2);
  
  total_tools ~ poisson(lambda);
}

generated quantities {

  // Posterior Predictive samples
  real sim_lambda[n];
  int sim_total_tools[n];
  
  // Loglikelihood
  real log_lik[n];
  
  for (i in 1:n) {
    sim_lambda[i] = exp(alpha[cid[i]] + beta[cid[i]] * P[i]);
    sim_total_tools[i] = poisson_rng(sim_lambda[i]);
    
    log_lik[i] = poisson_lpmf(total_tools[i] | lambda[i]);
  }
}
"

stan_poisson_m_1 <- stan_model(model_code = poisson_m1_code)

fit_stan_poisson_m_1 <- 
  sampling(stan_poisson_m_1, 
           iter = 1000, chains = 4, 
           data = c(compose_data(d %>% 
                                   select(total_tools, 
                                          contact, 
                                          P,
                                          cid = contact_id)),
                    N_contacts = length(unique(d$contact))),
           verbose = TRUE)

print(fit_stan_poisson_m_1, pars = c("alpha", "beta"))

the_pars <- c("alpha[1]", "alpha[2]", "beta[1]", "beta[2]")

mcmc_rank_hist(fit_stan_poisson_m_1, pars = )
mcmc_rank_overlay(fit_stan_poisson_m_1, pars = the_pars)
mcmc_pairs(fit_stan_poisson_m_1, pars = the_pars)
mcmc_intervals(fit_stan_poisson_m_1, pars = the_pars)
mcmc_areas(fit_stan_poisson_m_1, pars = the_pars)
bayesplot::mcmc_acf(fit_stan_poisson_m_1, pars = the_pars)

PSIS(fit_stan_poisson_m_1)
rethinking::WAIC(fit_stan_poisson_m_1)

loo(fit_stan_poisson_m_1)
loo::waic(extract(fit_stan_poisson_m_1, "log_lik")[[1]])


# PRIORS ------------------------------------------------------------------

# log(lambda_i) = alpha
# alpha ~ normal(0, 10) // FLAT

curve(dlnorm(x, 0, 10), from = 0, to = 100, n = 200,
      xlab = "Mean number of tools")

mean(exp(rnorm(1e4, 0, 10)))

curve(dlnorm(x, 3, 0.5), from = 0, to = 100, n = 200, col = "blue", add = TRUE)

mean(exp(rnorm(1e4, 3, 0.5)))

# Flat prior for beta
N <- 100
alpha <- rnorm(N, 3, 0.5)
beta <- rnorm(N, 0, 10)

plot(NULL, xlim = c(-2, 2), ylim = c(0, 100),
     ylab = "total_tools", xlab = "scaled log population",
     main = "beta ~ dnorm(0, 10")

for (i in 1:N) curve(exp(alpha[i] + beta[i] * x), add = TRUE, col=grau())

# Regularize
beta <- rnorm(N, 0, 0.2)

plot(NULL, xlim = c(-2, 2), ylim = c(0, 100),
     ylab = "total_tools", xlab = "scaled log population",
     main = "beta ~ dnorm(0, 0.2)")

for (i in 1:N) curve(exp(alpha[i] + beta[i] * x), add = TRUE, col=grau())

# Non-scaled log population
log_pop_seq <- seq(log(100), log(200000), length.out = N)

lambda <- sapply(log_pop_seq, function(x) exp(alpha + beta*x))

plot(NULL, xlim = range(log_pop_seq), ylim = c(0, 500),
     ylab = "total_tools", xlab = "log population",
     main = "beta ~ dnorm(0, 0.2)")
for (i in 1:N) lines(log_pop_seq, lambda[i, ], col = grau(), lwd = 1.5)

# Population
plot(NULL, xlim = range(exp(log_pop_seq)), ylim = c(0, 500),
     ylab = "total_tools", xlab = "population",
     main = "beta ~ dnorm(0, 0.2)")
for (i in 1:N) lines(exp(log_pop_seq), lambda[i, ], col = grau(), lwd = 1.5)


# POSTERIOR SAMPLIN -------------------------------------------------------

# Intercept-only model

get_prior(total_tools ~ -1 + factor(contact),
          data = d,
          family = poisson)

brm_poisson_m2_priors <- 
  c(prior(normal(3, 0.5), class = b))

make_stancode(total_tools ~ -1 + factor(contact) ,
              data = d,
              family = poisson,
              prior = brm_poisson_m2_priors)

fit_brm_poisson_m2 <- 
  brm(total_tools ~ -1 + factor(contact) ,
      data = d,
      family = poisson,
      prior = brm_poisson_m2_priors,
      iter = 1000, warmup = 500, chains = 4, cores = 4, seed = 2021)

fit_brm_poisson_m2_stan <- fit_brm_poisson_m2$fit

c_id <- sapply(1:4, rep, 500)
r_eff_m1 <- loo::relative_eff(extract(fit_stan_poisson_m_1, "log_lik")[[1]], 
                              chain_id = c_id)
r_eff_m2 <- loo::relative_eff(log_lik(fit_brm_poisson_m2), 
                              chain_id = c_id)

loo_m1 <- loo::loo(extract(fit_stan_poisson_m_1, "log_lik")[[1]], 
                   r_eff = r_eff_m1)
loo_m2 <- loo::loo(log_lik(fit_brm_poisson_m2), 
                   r_eff = r_eff_m1)

loo_comp <- loo::loo_compare(loo_m1, loo_m2)

# PPS
post_lambda <- post$sim_lambda %>% 
  apply(2, tidybayes::mean_hdci, .width = 0.89) %>% 
  bind_rows()

post_total_tools <- post$sim_total_tools %>% 
  apply(2, tidybayes::mean_hdci, .width = 0.89) %>% 
  bind_rows()

plot(NULL, xlim = range(d$population), ylim = c(0,100))
lines(total_tools ~ population, d, type = "p")

lines(d$population, post_lambda$y)
lines(d$population, post_lambda$ymin, lty = 2)
lines(d$population, post_lambda$ymax, lty = 2)

lines(d$population, post_total_tools$y, col = "blue")
lines(d$population, post_total_tools$ymin, col = "blue", lty = 2)
lines(d$population, post_total_tools$ymax, col = "blue", lty = 2)

###
plot(total_tools ~ P, d, type = "p")

lines(d$P, post_lambda$y)
lines(d$P, post_lambda$ymin, lty = 2)
lines(d$P, post_lambda$ymax, lty = 2)

lines(d$P, post_total_tools$y, col = "blue")
lines(d$P, post_total_tools$ymin, col = "blue", lty = 2)
lines(d$P, post_total_tools$ymax, col = "blue", lty = 2)

## BIEN hecho
post <- extract(fit_stan_poisson_m_1, pars = c("alpha", "beta"))
dim(post$alpha)
dim(post$beta)

idx_samples <- sample(1:2000, 1e4, replace = TRUE)

post_alpha <- post$alpha[idx_samples, ]
post_beta <- post$beta[idx_samples, ]

post <- extract(fit_stan_poisson_m_1, pars = c("sim_lambda", "sim_total_tools"))

P_seq <- seq(range(d$P)[1], range(d$P)[2], length.out = 1000)

# cid = 1 (LOW contact)
post_lambda <- sapply(P_seq,
                      function(p) exp(post_alpha[,1] + 
                                        post_beta[, 1] * p))
post_lambda_cis <- post_lambda %>% 
  apply(2, tidybayes::mean_qi, .width = 0.89) %>% 
  bind_rows()

plot(total_tools ~ P, d, type = "p")

lines(P_seq, post_lambda_cis$y)
lines(P_seq, post_lambda_cis$ymin, lty = 2)
lines(P_seq, post_lambda_cis$ymax, lty = 2)

# cid = 2 (HIGH contact)
post_lambda <- sapply(P_seq,
                      function(p) exp(post_alpha[, 2] + 
                                        post_beta[, 2] * p))
post_lambda_cis <- post_lambda %>% 
  apply(2, tidybayes::mean_qi, .width = 0.89) %>% 
  bind_rows()

lines(P_seq, post_lambda_cis$y, col = "blue")
lines(P_seq, post_lambda_cis$ymin, col = "blue", lty = 2)
lines(P_seq, post_lambda_cis$ymax, col = "blue", lty = 2)

