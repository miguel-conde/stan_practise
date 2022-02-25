library(tidyverse)
library(rethinking)
library(brms)
library(tidybayes)
library(bayesplot)
library(rstan)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

source(here::here("rethinking", "help_functions.R"), encoding = "utf-8")

# Introducció a los modelos multinivel.
# 
# Estrategia: partial pooling mediante varying intercepts
#


# 1 - Modelo RENACUAJOS ---------------------------------------------------


# Data --------------------------------------------------------------------


data(reedfrogs)

d <- reedfrogs %>% as_tibble()

d

help(reedfrogs)

glimpse(d)

Hmisc::describe(d)

psych::describe(d)

psych::describeBy(d, d$size)

psych::describeBy(d, d$pred)

# Nos interesa 'surv', numero de renacuajos supervivientes de los 'density' 
# iniciales

d %>% group_by(density, pred, size) %>% summarise(n = n())

# Hay 4 "individuos" (tanks) en cada uno de los 12 grupos (clusters) existentes


# 2 - Modelo CHIMPANCÉS ---------------------------------------------------


# Data --------------------------------------------------------------------



data(chimpanzees)

d <- chimpanzees %>% as_tibble()

d

help(chimpanzees)

glimpse(d)

Hmisc::describe(d)

psych::describe(d)

psych::describeBy(d, d$actor)

psych::describeBy(d, d$block)

d %>% group_by(actor) %>% summarise(n = n())
d %>% group_by(block) %>% summarise(n = n())
d %>% group_by(actor, block) %>% summarise(n = n())
d %>% group_by(block, actor) %>% summarise(n = n())

d <- d %>% mutate(treatment = as.integer(1 + prosoc_left + 2 * condition))

d %>% group_by(treatment) %>% summarise(n = n())
d %>% group_by(block, actor, treatment) %>% summarise(n = n())

# 2 tipos de cluster ------------------------------------------------------

stan_m1_code <- "
data {
  int<lower = 1> N_actors;
  int<lower = 1> N_blks;
  int<lower = 1> N_treatments;
  int<lower = 1> n;
  int<lower = 1, upper = N_actors> actor[n];
  int<lower = 1, upper = N_blks> blk[n];
  int<lower = 1, upper = N_treatments> treatment[n];
  int<lower = 0, upper = 1> pulled_left[n];
}

parameters {
  real alpha[N_actors];
  real gamma[N_blks];
  real beta[N_treatments];
  
  real bar_alpha;
  real<lower = 0> sigma_alpha;
  real<lower = 0> sigma_gamma;

}

transformed parameters {
  vector[n] p;
  
  for (i in 1:n) {
    p[i] =  alpha[actor[i]] + gamma[blk[i]] + beta[treatment[i]];
    p[i] =  inv_logit(p[i]);
  }

}

model {  

  // PRIORS
  alpha ~ normal(bar_alpha, sigma_alpha);
  gamma ~ normal(0, sigma_gamma);
  beta ~ normal(0, 0.5);
  
  // HYPERPRIORS
  bar_alpha ~ normal(0, 1.5);
  sigma_alpha ~ exponential(1);
  sigma_gamma ~ exponential(1);
  
  pulled_left ~ binomial(1, p);
}

generated quantities {

}
"

stan_m1 <- stan_model(model_code = stan_m1_code)

fit_stan_m1 <- sampling(stan_m1, 
                        iter = 2000, chains = 1, 
                        control = list(adapt_delta = 0.99),
                        data = c(compose_data(d %>% 
                                                select(actor,
                                                       blk = block,
                                                       treatment,
                                                       pulled_left)), 
                                 N_actors = length(unique(d$actor)),
                                 N_blks = length(unique(d$block)),
                                 N_treatments = length(unique(d$treatment))),
                        verbose = TRUE)

the_pars <- c("alpha", "gamma", "beta", 
              "bar_alpha", "sigma_alpha", "sigma_gamma")
print(fit_stan_m1, pars = the_pars)
plot(fit_stan_m1, pars = the_pars)
# rethinking::pairs(fit_stan_m1, pars = the_pars)


# NON CENTERED REFACTORING ------------------------------------------------

## Para solucionear DIVERGENCIAS

stan_m2_code <- "
data {
  int<lower = 1> N_actors;
  int<lower = 1> N_blks;
  int<lower = 1> N_treatments;
  int<lower = 1> n;
  int<lower = 1, upper = N_actors> actor[n];
  int<lower = 1, upper = N_blks> blk[n];
  int<lower = 1, upper = N_treatments> treatment[n];
  int<lower = 0, upper = 1> pulled_left[n];
}

parameters {
  vector[N_actors] x;
  vector[N_blks] z;
  vector[N_treatments] beta;
  
  real bar_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;

}

transformed parameters { 
  vector[n] p;
  
  for (i in 1:n) {
    p[i] =  bar_alpha + sigma_alpha * x[actor[i]] + sigma_gamma * z[blk[i]] + beta[treatment[i]];
    p[i] =  inv_logit(p[i]);
  }

}

model {  
  
  // HYPERPRIORS
  bar_alpha ~ normal(0, 1.5);
  sigma_alpha ~ exponential(1);
  sigma_gamma ~ exponential(1);

  // PRIORS
  x ~ normal(0, 1);
  z ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  
  pulled_left ~ binomial(1, p);
}

generated quantities {

}
"

stan_m2 <- stan_model(model_code = stan_m2_code)

fit_stan_m2 <- sampling(stan_m2, 
                        iter = 1000, chains = 4, 
                        data = c(compose_data(d %>% 
                                                select(actor,
                                                       blk = block,
                                                       treatment,
                                                       pulled_left)), 
                                 N_actors = length(unique(d$actor)),
                                 N_blks = length(unique(d$block)),
                                 N_treatments = length(unique(d$treatment))),
                        verbose = TRUE)

print(fit_stan_m2, pars = the_pars)
plot(fit_stan_m2, pars = the_pars)

## Lo mismo pero conservando los coeficientes originales

stan_m3_code <- "
data {
  int<lower = 1> N_actors;
  int<lower = 1> N_blks;
  int<lower = 1> N_treatments;
  int<lower = 1> n;
  int<lower = 1, upper = N_actors> actor[n];
  int<lower = 1, upper = N_blks> blk[n];
  int<lower = 1, upper = N_treatments> treatment[n];
  int<lower = 0, upper = 1> pulled_left[n];
}

parameters {
  vector[N_actors] x;
  vector[N_blks] z;
  vector[N_treatments] beta;
  
  real bar_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;

}

transformed parameters {
  vector[N_actors] alpha;
  vector[N_blks] gamma;
  vector[n] p;
  
  for (i in 1:n) {
  
    alpha[actor[i]] = bar_alpha + sigma_alpha * x[actor[i]];
    gamma[blk[i]] = sigma_gamma * z[blk[i]];
    
    p[i] =  alpha[actor[i]] + gamma[blk[i]] + beta[treatment[i]];
    p[i] =  inv_logit(p[i]);
  }

}

model {  
  
  // HYPERPRIORS
  bar_alpha ~ normal(0, 1.5);
  sigma_alpha ~ exponential(1);
  sigma_gamma ~ exponential(1);

  // PRIORS
  x ~ normal(0, 1);
  z ~ normal(0, 1);
  beta ~ normal(0, 0.5);
  
  pulled_left ~ binomial(1, p);
}

generated quantities {

  vector[n] log_lik;
  
  for (i in 1:n) {
  
    log_lik[i] = binomial_lpmf(pulled_left | 1, p[i]);
  }

}
"

stan_m3 <- stan_model(model_code = stan_m3_code)

fit_stan_m3 <- sampling(stan_m3, 
                        iter = 4000, chains = 1, 
                        data = c(compose_data(d %>% 
                                                select(actor,
                                                       blk = block,
                                                       treatment,
                                                       pulled_left)), 
                                 N_actors = length(unique(d$actor)),
                                 N_blks = length(unique(d$block)),
                                 N_treatments = length(unique(d$treatment))),
                        verbose = TRUE)

print(fit_stan_m3, pars = the_pars)
plot(fit_stan_m3, pars = the_pars)

traceplot(fit_stan_m3, pars = the_pars)


# Con Ulam ----------------------------------------------------------------


dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  blk = d$block,
  treatment = as.integer(d$treatment))


set.seed(13)
m13.4nc <- ulam(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a_bar + z[actor] * sigma_a + # actorintercepts
      x[blk] * sigma_g +                     # blockintercepts
      b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    z[actor] ~ dnorm(0, 1),
    x[blk] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    gq> vector[actor]:a <<- a_bar + z * sigma_a,
    gq> vector[blk]:g <<- x * sigma_g
  ), data = dat_list, chains = 4, cores = 4)


# Con BRMS ----------------------------------------------------------------


detach(package:rethinking, unload = T)

bf <- pulled_left | trials(1) ~ 1 + (1 | actor) + (1 | block) + (1 | treatment)

get_prior(bf,
          data = d,
          family = binomial)

priors_m3 <- 
  c(prior(normal(0, 1.5), class = Intercept),
    prior(exponential(1), class = sd, coef = Intercept, group = actor),
    prior(exponential(1), class = sd, coef = Intercept, group = block),
    prior(normal(0.5, 1e-16), class = sd, coef = Intercept, group = treatment))


make_stancode(formula = bf,
              data = d,
              family = binomial,
              prior = priors_m3)

fit_brm_m3 <- brm(formula = bf,
                  data = d,
                  family = binomial,
                  prior = priors_m3,
                  iter = 1000, warmup = 500, chains = 1, cores = 4, seed = 2021)

fit_brm_binom_m1_stan <- fit_brm_m3$fit

fit_brm_m3
pairs(fit_brm_m3)
plot(fit_brm_m3)

coef(fit_brm_m3)
fixef(fit_brm_m3)
ranef(fit_brm_m3)


# 3 - VARYING EFFECTS -----------------------------------------------------

library(rethinking)

a <- 3         # average morning wait time
b <- (-1)      # average difference after noon wait time
sigma_a <- 1   # std dev in intercepts
sigma_b <- 0.5 # std dev in slopes
rho <- (-0.7)  # correlation between intercepts and slopes

Mu <- c(a, b)

sigmas <- c(a, b)
Sigmas <- diag(sigmas)

Rho <- matrix(c(1, rho, rho, 1), nrow = 2)

Sigma <- Sigmas %*% Rho %*% Sigmas

N_cafes <- 20
set.seed(5)
vary_effects <- MASS::mvrnorm(N_cafes, Mu, Sigma)
dimnames(vary_effects) <- list(n = NULL, var = c("Intercept", "slope"))

plot(vary_effects[, "Intercept"], vary_effects[, "slope"], col = rangi2,
      xlab= "intercepts (a_cafe)", ylab="slopes(b_cafe)")

# overlay population distribution

library(ellipse)

for (l in c(0.1, 0.3, 0.5, 0.8, 0.99))
  lines(ellipse(Sigma, centre = Mu, level = l), col = col.alpha("black", 0.2))

set.seed(22) 

N_visits <- 10

afternoon <- rep(0:1, N_visits * N_cafes / 2)
cafe_id <- rep(1:N_cafes, each = N_visits)

mu <- vary_effects[cafe_id, "Intercept"] + vary_effects[cafe_id, "slope"] * afternoon
sigma <- 0.5 # std dev within cafes

wait <- rnorm(N_visits * N_cafes, mu, sigma)

d <- tibble(cafe = cafe_id, afternoon = afternoon, wait = wait)

##
R <- rlkjcorr(1e4, K = 2, eta = 2)
dens(R[, 1, 2], xlab = "correlation")

a_cafe <- vary_effects[, "Intercept"]
b_cafe <- vary_effects[, "slope"]

set.seed(867530)
m14.1 <- ulam(
  alist(
    wait ~ normal(mu, sigma),
    mu <- a_cafe[cafe] + b_cafe[cafe] * afternoon,
    c(a_cafe, b_cafe)[cafe] ~ multi_normal(c(a, b), Rho, sigma_cafe),
    a ~ normal(5,2),
    b ~ normal(-1,0.5),
    sigma_cafe ~ exponential(1),
    sigma ~ exponential(1),
    Rho ~ lkj_corr(2)
  ), data = d, chains = 4, cores = 4)

rethinking::stancode(m14.1)


# 13H1 --------------------------------------------------------------------

data(bangladesh)

d <- bangladesh %>% as_tibble()

glimpse(d)
summary(d)

d <- d %>% mutate(district = as.integer(factor(district))) %>% janitor::clean_names()

stan_code_13H1_pooled <- "
data {
  int<lower = 1> n;
  
  int<lower = 1> district[n];
  int<lower = 0, upper = 1> use_contraception[n];
}

parameters {
  real alpha;
}

transformed parameters {

  vector[n] p;
  
  for (i in 1:n) {
    p[i] = alpha;
    p[i] = inv_logit(p[i]);
  }

}

model {

  alpha ~ normal(0, 1.5);

  use_contraception ~ binomial(1, p);

}

generated quantities {

  vector[n] log_lik;
  
  for (i in 1:n) {
    log_lik[i] = binomial_lpmf(use_contraception | 1, p[i]);
  }

}

"

stan_m1_pooled <- stan_model(model_code = stan_code_13H1_pooled)

fit_stan_m1_pooled <- sampling(stan_m1_pooled, 
                               iter = 2000, chains = 1, 
                               # control = list(adapt_delta = 0.99),
                               data = c(compose_data(d %>% 
                                                       select(district, 
                                                              use_contraception))),
                               verbose = TRUE)

the_pars <- c("alpha")
print(fit_stan_m1_pooled, pars = the_pars)
plot(fit_stan_m1_pooled, pars = the_pars)

inv_logit(extract(fit_stan_m1_pooled, pars = "alpha")[[1]]) %>% mean()
inv_logit(extract(fit_stan_m1_pooled, pars = "alpha")[[1]]) %>% tidybayes::median_qi()
d$use_contraception %>% mean()

# Only fixed-effects

stan_code_13H1_non_pooled <- "
data {
  int<lower = 1> n;
  int<lower = 1> N_districts;
  
  int<lower = 1> district[n];
  int<lower = 0, upper = 1> use_contraception[n];
}

parameters {
  vector[N_districts] alpha;
}

transformed parameters {

  vector[n] p;
  
  for (i in 1:n) {
    p[i] = alpha[district[i]];
    p[i] = inv_logit(p[i]);
  }

}

model {

  alpha ~ normal(0, 1.5);

  use_contraception ~ binomial(1, p);

}

generated quantities {

  vector[n] log_lik;
  
  for (i in 1:n) {
    log_lik[i] = binomial_lpmf(use_contraception | 1, p[i]);
  }

}

"

stan_m1_non_pooled <- stan_model(model_code = stan_code_13H1_non_pooled)

fit_stan_m1_non_pooled <- sampling(stan_m1_non_pooled, 
                                   iter = 2000, chains = 1, 
                                   # control = list(adapt_delta = 0.99),
                                   data = c(compose_data(d %>% 
                                                           select(district, 
                                                                  use_contraception)),
                                            N_districts = d$district %>% unique() %>% length()),
                                   verbose = TRUE)

the_pars <- c("alpha")
print(fit_stan_m1_non_pooled, pars = the_pars)
plot(fit_stan_m1_non_pooled, pars = the_pars)

inv_logit(extract(fit_stan_m1_non_pooled, pars = "alpha")[[1]]) %>% 
  apply(2, mean)
inv_logit(extract(fit_stan_m1_non_pooled, 
                  pars = "alpha")[[1]]) %>% apply(2, tidybayes::median_qi) %>% 
  bind_rows()
d %>% group_by(district) %>% summarise(p = mean(use_contraception))

# Varying effects  -only intercetps

stan_code_13H1_semi_pooled <- "
data {
  int<lower = 1> n;
  int<lower = 1> N_districts;
  
  int<lower = 1> district[n];
  int<lower = 0, upper = 1> use_contraception[n];
}

parameters {
  real alpha_bar;
  real<lower = 0> sigma;
  vector[N_districts] alpha;
}

transformed parameters {

  vector[n] p;
  
  for (i in 1:n) {
    p[i] = alpha[district[i]];
    p[i] = inv_logit(p[i]);
  }

}

model {

  alpha_bar ~ normal(0, 1.5);
  sigma ~ exponential(1);
  
  alpha ~ normal(alpha_bar, sigma);

  use_contraception ~ binomial(1, p);

}

generated quantities {

  vector[n] log_lik;
  
  for (i in 1:n) {
    log_lik[i] = binomial_lpmf(use_contraception | 1, p[i]);
  }

}

"

stan_m1_semi_pooled <- stan_model(model_code = stan_code_13H1_semi_pooled)

fit_stan_m1_semi_pooled <- sampling(stan_m1_semi_pooled, 
                                   iter = 2000, chains = 1, 
                                   # control = list(adapt_delta = 0.99),
                                   data = c(compose_data(d %>% 
                                                           select(district, 
                                                                  use_contraception)),
                                            N_districts = d$district %>% unique() %>% length()),
                                   verbose = TRUE)

the_pars <- c("alpha", "alpha_bar", "sigma")
print(fit_stan_m1_semi_pooled, pars = the_pars)
plot(fit_stan_m1_semi_pooled, pars = the_pars)

inv_logit(extract(fit_stan_m1_semi_pooled, 
                  pars = "alpha")[[1]]) %>% apply(2, tidybayes::median_qi) %>% 
  bind_rows()

inv_logit(extract(fit_stan_m1_semi_pooled, 
                  pars = "alpha_bar")[[1]]) %>% tidybayes::median_qi()

## Plot

pooled_fits <- extract(fit_stan_m1_pooled, pars = "alpha")[[1]] %>% 
  inv_logit() %>% tidybayes::median_qi(.width = 0.89)

non_pooled_fits <- extract(fit_stan_m1_non_pooled, pars = "alpha")[[1]] %>% 
  inv_logit() %>% 
  apply(2, tidybayes::median_qi, .width = 0.89) %>% 
  bind_rows()

semi_pooled_fits <- extract(fit_stan_m1_semi_pooled, pars = "alpha")[[1]] %>% 
  inv_logit() %>% 
  apply(2, tidybayes::median_qi, .width = 0.89) %>% 
  bind_rows()

probe_plot <- d %>% group_by(district) %>% 
  summarise(n = n(), use_contraception = mean(use_contraception)) %>% 
  arrange(n)

plot(use_contraception ~ district, probe_plot, pch = 0)
abline(h = pooled_fits$y, lty = 2, col = "black", pch = 1)
lines(probe_plot$district, non_pooled_fits$y, col = "blue", type = "p", pch = 2)
lines(probe_plot$district, semi_pooled_fits$y, col = "red", type = "p", pch = 3)
legend("topright",
       legend = c("Raw", "Pooled", "Non pooled", "Semi-pooled"),
       pch = c(0, NA, 2, 3), 
       lty = c(0, 2, 0, 0),
       col = c("black", "black", "blue", "red"),
       bty = "n")

plot(use_contraception ~ n, probe_plot, pch = 0)
abline(h = pooled_fits$y, lty = 2, col = "black", pch = 1)
lines(probe_plot$n, non_pooled_fits$y, col = "blue", type = "p", pch = 2)
lines(probe_plot$n, semi_pooled_fits$y, col = "red", type = "p", pch = 3)
legend("topright",
       legend = c("Raw", "Pooled", "Non pooled", "Semi-pooled"),
       pch = c(0, NA, 2, 3), 
       lty = c(0, 2, 0, 0),
       col = c("black", "black", "blue", "red"),
       bty = "n")

waic_comp <- loo::loo_compare(loo::waic(loo::extract_log_lik(fit_stan_m1_pooled)),
                              loo::waic(loo::extract_log_lik(fit_stan_m1_non_pooled)),
                              loo::waic(loo::extract_log_lik(fit_stan_m1_semi_pooled)))

print(waic_comp, simplify = FALSE)

psis_comp <- loo::loo_compare(loo::loo(loo::extract_log_lik(fit_stan_m1_pooled)),
                           loo::loo(loo::extract_log_lik(fit_stan_m1_non_pooled)),
                           loo::loo(loo::extract_log_lik(fit_stan_m1_semi_pooled)))

print(psis_comp, simplify = FALSE)

tibble(value             = c("difference", "se"),
       elpd              = m_comp[3, 1:2],
       conversion_factor = c(-2, 2)) %>% 
  mutate(waic            = elpd * conversion_factor)

rethinking::compare(fit_stan_m1_pooled,
                    fit_stan_m1_non_pooled, 
                    fit_stan_m1_semi_pooled)
rethinking::compare(fit_stan_m1_pooled,
                    fit_stan_m1_non_pooled, 
                    fit_stan_m1_semi_pooled,
                    func = "PSIS")

compute_lppd(extract(fit_stan_m1_semi_pooled, "log_lik")[[1]])
compute_pWAIC(extract(fit_stan_m1_semi_pooled, "log_lik")[[1]])
compute_WAIC(extract(fit_stan_m1_semi_pooled, "log_lik")[[1]]
             )
