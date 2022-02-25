library(tidyverse)
library(rethinking)
library(brms)
library(tidybayes)
library(bayesplot)
library(rstan)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


# DATA --------------------------------------------------------------------


data(chimpanzees)

d <- chimpanzees %>% as_tibble()

# Response: pulled_left - El lado de la palanca que movió el chimpancé
#
# Predictors:
#              prosoc_left: en qué lado estaba la opción prosocial
#              condition: si había o no otro chimpancé al otro lado de la mesa

# Variable índice
#
# 1 - prosoc_left = 0 - condition = 0 2 food right / no partner
# 2 - prosoc_left = 1 - condition = 0 2 food right / no partner
# 3 - prosoc_left = 0 - condition = 1 2 food left / partner present
# 4 - prosoc_left = 1 - condition = 1 2 food left / partner present
#

d <- d %>% mutate(treatment = 1 + prosoc_left + 2 * condition)

xtabs( ~ treatment + prosoc_left + condition, d)


# PRIORS ------------------------------------------------------------------

## Prior predictive sampling

# Flat prior for alpha (actor)

prior_alpha <- rnorm(1e5, 0, 10)
prior_p <- inv_logit(prior_alpha)
prior_pulled_left <- rbinom(1e5, 1, prior_p)

plot(density(prior_p))

# Regularize
prior_alpha <- rnorm(1e5, 0, 1.5)
prior_p <- inv_logit(prior_alpha)
prior_pulled_left <- rbinom(1e5, 1, prior_p)

lines(density(prior_p), col = "blue")

# Let's see beta (treatment)
prior_beta <- sapply(1:4, function(x) rnorm(1e5, 0, 10))

prior_p <- sapply(1:4, function(x) inv_logit(prior_alpha + 
                    prior_beta[, x]))

plot(density(prior_p))
plot(density(abs(prior_p[, 1] - prior_p[, 2])))

# Regularize
prior_beta <- sapply(1:4, function(x) rnorm(1e5, 0, 0.5))

prior_p <- sapply(1:4, function(x) inv_logit(prior_alpha + 
                                               prior_beta[, x]))

lines(density(abs(prior_p[, 1] - prior_p[, 2])), col = "blue")

mean(abs(prior_p[, 1] - prior_p[, 2]))

# Ahora el modelo "cree" a priori en que hay poca diferencia entre tratamientos


# FIT ---------------------------------------------------------------------

# Stan --------------------------------------------------------------------

binom_m_1_code <- "
data {
  int<lower = 1> N_actors;
  int<lower = 1> N_treatments;
  int<lower = 1> n;
  int<lower = 1, upper = N_treatments> treatment[n]; // Predictor
  int<lower = 1, upper = N_actors>     actor[n];     // Predictor
  int<lower = 0, upper = 1> pulled_left[n];          // Response
  
}

parameters {
  real alpha[N_actors];
  real beta[N_treatments];
}

transformed parameters {
  vector[n] p;
  
  for (i in 1:n) {
    p[i] =  alpha[actor[i]] + beta[treatment[i]];
    p[i] =  inv_logit(p[i]);
  }
}

model {

  pulled_left ~ binomial(1, p); // L ~ bernouilli(p);
  
  alpha ~ normal(0, 1.5);
  beta ~ normal(0, 0.5);
  
  pulled_left ~ binomial(1, p);
}

generated quantities {
  vector[n] log_lik;
  vector[n] sim_p;
  int sim_pulled_left[n];
  
  for (i in 1:n) {
  
    sim_p[i] =  alpha[actor[i]] + beta[treatment[i]];
    sim_p[i] =  inv_logit(p[i]);
    
    sim_pulled_left[i] = binomial_rng(1, sim_p[i]);
    
    log_lik[i] = binomial_lpmf(pulled_left[i] | 1, p[i]);
  }
}
"

stan_binom_m_1 <- stan_model(model_code = binom_m_1_code)

fit_stan_binom_m_1 <- sampling(stan_binom_m_1, 
                          iter = 1000, chains = 1, 
                          data = c(compose_data(d %>% 
                                                  select(treatment, 
                                                         actor, 
                                                         pulled_left)),
                                   N_treatments = length(unique(d$treatment)), 
                                   N_actors = length(unique(d$actor))),
                          verbose = TRUE)

print(fit_stan_binom_m_1, pars = c("alpha", "beta"))


# Ulam --------------------------------------------------------------------

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment))

fit_ulam_binom_m_1 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data = dat_list, chains = 4, log_lik = TRUE, cmdstan = TRUE)

precis(fit_ulam_binom_m_1, depth = 2)


# brms --------------------------------------------------------------------

get_prior(pulled_left | trials(1) ~ -1 + factor(actor) + factor(treatment),
          data = d,
          family = binomial)

brm_binom_m1_priors <- 
  c(prior(normal(0, 1.5), class = b))

make_stancode(pulled_left | trials(1) ~ -1 + factor(actor) + factor(treatment),
              data = d,
              family = binomial,
              prior = brm_binom_m1_priors)

fit_brm_binom_m1 <- brm(pulled_left | trials(1) ~ -1 + factor(actor) + factor(treatment),
                        data = d,
                        family = binomial,
                        prior = brm_binom_m1_priors,
                        iter = 1000, warmup = 500, chains = 4, cores = 4, seed = 2021)

fit_brm_binom_m1_stan <- fit_brm_binom_m1$fit

fit_brm_binom_m1
pairs(fit_brm_binom_m1)
plot(fit_brm_binom_m1)

coef(fit_brm_binom_m1)

###

get_prior(pulled_left | trials(1) ~ -1 + (1 | actor) + (1 | treatment),
          data = d,
          family = binomial)

brm_binom_m1_priors <- 
  c(prior(normal(0, 1.5), class = sd, coef = Intercept, group = actor),
    prior(normal(0, 0.5), class = sd, coef = Intercept, group = treatment))

make_stancode(pulled_left | trials(1) ~ -1 + (1 | actor) + (1 | treatment),
              data = d,
              family = binomial,
              prior = brm_binom_m1_priors)

fit_brm_binom_m1 <- brm(pulled_left | trials(1) ~ -1 + (1 | actor) + (1 | treatment),
                        data = d,
                        family = binomial,
                        prior = brm_binom_m1_priors,
                        iter = 1000, warmup = 500, chains = 4, cores = 4, seed = 2021)

fit_brm_binom_m1_stan <- fit_brm_binom_m1$fit

fit_brm_binom_m1
pairs(fit_brm_binom_m1)
plot(fit_brm_binom_m1)

coef(fit_brm_binom_m1)


# POSTERIOR ---------------------------------------------------------------

## PPS
range_idx <- 1 : (fit_stan_binom_m_1@stan_args[[1]]$iter - 
  fit_stan_binom_m_1@stan_args[[1]]$warmup) 

idx_sample <- sample(range_idx, 1e5, replace = TRUE)
post_alpha <- extract(fit_stan_binom_m_1, "alpha")[[1]][idx_sample, ] %>% 
  apply(2, inv_logit)

boxplot(post_alpha)
post_alpha %>% apply(2, tidybayes::median_qi, .width = .89) %>% bind_rows()
post_alpha %>% apply(2, tidybayes::median_hdci, .width = .89) %>% bind_rows()

# > .5 => la predicción sería prosoc_left
# < .5 => la predicción sería prosoc_right

post_beta <- extract(fit_stan_binom_m_1, "beta")[[1]][idx_sample, ] %>% 
  apply(2, inv_logit)

boxplot(post_beta)
post_beta %>% apply(2, tidybayes::median_qi) %>% bind_rows()
post_beta %>% apply(2, tidybayes::median_hdci) %>% bind_rows()

# Contrastes (scala de probabilidades):
boxplot(post_beta[,1] - post_beta[,3],
        post_beta[,2] - post_beta[,4])

# Contrastes (escala de log odds):
post_beta_log_odds <- extract(fit_stan_binom_m_1, "beta")[[1]][idx_sample, ]

boxplot(post_beta_log_odds[,1] - post_beta_log_odds[,3],
        post_beta_log_odds[,2] - post_beta_log_odds[,4])

## PP Check

d %>% 
  group_by(actor, treatment) %>% 
  summarise(pulled_left = mean(pulled_left)) %>% 
  spread(treatment, pulled_left)

post_alpha <- extract(fit_stan_binom_m_1, "alpha")[[1]][idx_sample, ]
post_beta  <- extract(fit_stan_binom_m_1, "beta")[[1]][idx_sample, ]

post_p <- sapply(1:nrow(d), function(i) {
  inv_logit(post_alpha[,d$actor[i]] + post_beta[,d$treatment[i]])
}) 

post_pulled_left <- post_p %>% 
  apply(1, function(p) rbinom(n = length(p), size = 1, p))

# Con PPS de p
p_hdpi <- post_p %>% apply(2, tidybayes::median_hdci, .width = 0.89) %>% 
  bind_rows()
p_hdpi

pred_d <- d %>% select(actor, treatment) %>% 
  bind_cols(as_tibble(p_hdpi)) 

pred_d %>% 
  group_by(actor, treatment) %>% 
  summarise(y = mean(y)) %>% 
  spread(treatment, y)

# Con PPS pulled left
d %>% select(actor, treatment) %>%
  mutate(mean_pulled_left = post_pulled_left %>% apply(1, mean)) %>% 
  group_by(actor, treatment) %>% 
  summarise(pulled_left = mean(mean_pulled_left)) %>% 
  spread(treatment, pulled_left)

