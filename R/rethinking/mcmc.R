library(tidyverse)

library(rethinking)

library(tidybayes)

library(dagitty)

library(brms)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


# 9H1 ---------------------------------------------------------------------

# Run the model below and then inspect the posterior distribution and explain 
# what it is accomplishing.

mp <- ulam(
  alist(
    a ~ dnorm(0, 1),
    b ~ dcauchy(0, 1)
  ), data = list(y = 1), chains = 1)

# Compare the samples for the parameters a and b. Can you explain the 
# different traceplots? If you are unfamiliar with the Cauchy distribution, you
# should look it up. The key feature to attend to is that it has no expected 
# value. Can you connect this fact to the traceplot?

print(mp)
precis(mp)

traceplot(mp)


# 9H2 ---------------------------------------------------------------------

# Recall the divorce rate example from Chapter 5. Repeat that analysis, using 
# ulam this time, fitting models m5.1, m5.2, and m5.3. Use compare to compare
# the models on the basis of WAIC or PSIS.To use WAIC or PSIS with ulam, you 
# need add the argument log_log = TRUE. Explain the model comparison results.

data(WaffleDivorce)

d <- WaffleDivorce

d <- d %>% mutate(D = standardize(Divorce),
                  M = standardize(Marriage),
                  A = standardize(MedianAgeMarriage)) %>% 
  select(D, M, A)

m5.1 <- ulam(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), 
  data = d,
  chains = 4,
  log_lik = TRUE)

print(m5.1)
precis(m5.1)
pairs(m5.1)
pairs(m5.1@stanfit, pars = c("a", "bA", "sigma"))
traceplot(m5.1, pars = c("a", "bA", "sigma"))
trankplot(m5.1, pars = c("a", "bA", "sigma"))

m5.2 <- ulam(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), 
  data = d,
  chains = 4,
  log_lik = TRUE)

print(m5.2)
precis(m5.2)
pairs(m5.2)
pairs(m5.2@stanfit, pars = c("a", "bM", "sigma"))
traceplot(m5.2, pars = c("a", "bM", "sigma"))
trankplot(m5.2, pars = c("a", "bM", "sigma"))

m5.3 <- ulam(
  alist(
    D ~ dnorm(mu, sigma),
    mu <-a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), 
  data = d,
  chains = 4,
  log_lik = TRUE)

print(m5.3)
precis(m5.3)
pairs(m5.3)
pairs(m5.3@stanfit, pars = c("a", "bA", "bM", "sigma"))
traceplot(m5.3, pars = c("a", "bA", "bM", "sigma"))
trankplot(m5.3, pars = c("a", "bA", "bM", "sigma"))

compare(m5.1, m5.2, m5.3)
rethinking::WAIC(m5.1@stanfit)
rethinking::WAIC(m5.1@stanfit, pointwise = TRUE)

compare(m5.1, m5.2, m5.3, func = PSIS)

rethinking::PSIS(m5.1@stanfit)
rethinking::PSIS(m5.1@stanfit, pointwise = TRUE)

loo::loo(loo::extract_log_lik(m5.1@stanfit)) ## PSIS
res_psis <- loo::psis(loo::extract_log_lik(m5.1@stanfit))
res_psis

loo::waic(loo::extract_log_lik(m5.1@stanfit))

loo::loo_compare(loo::loo(loo::extract_log_lik(m5.1@stanfit)),
                 loo::loo(loo::extract_log_lik(m5.2@stanfit)),
                 loo::loo(loo::extract_log_lik(m5.3@stanfit)))

### brms

get_prior(D ~ A + 1,
          data = d,
          family = gaussian)

brm_m5.1_priors <- c(prior(normal(0, 0.2), class = Intercept),
                     prior(normal(0, 0.5), class = b, coef = A),
                     prior(exponential(1), class = sigma))

make_stancode(D ~ A + 1,
              data = d,
              family = gaussian,
              prior = brm_m5.1_priors)

brm_m5.1 <- brm(D ~ A + 1,
                data = d,
                family = gaussian,
                prior = brm_m5.1_priors,
                iter = 1000, warmup = 500, chains = 4, cores = 4, seed = 2021)
  
brm_m5.1
pairs(brm_m5.1)
plot(brm_m5.1)

### BAYESPLOT

library(bayesplot)
library(hrbrthemes)

post <- posterior_samples(brm_m5.1, add_chain = T)

mcmc_trace(post)

mcmc_trace(post[, c(1:3, 5)],  # we need to include column 7 because it contains the chain info 
           facet_args = list(nrow = 2, ncol = 2), 
           size = .15) +
  labs(title = "My custom trace plots") +
  scale_color_ipsum() +
  theme_ipsum() +
  theme(legend.position = c(.95, .2))

waic(brm_m5.1)
loo(brm_m5.1)

mcmc_acf(post, 
         pars = c("b_Intercept", "b_A", "sigma"),
         lags = 5) +
  scale_color_ipsum() +
  theme_ipsum()

###

brm_m5.1 <- add_criterion(brm_m5.1, c("waic", "loo", "bayes_R2", "loo_R2"))

brm_m5.1$criteria$waic
brm_m5.1$criteria$loo
brm_m5.1$criteria$bayes_R2 %>% tidybayes::median_qi()         # Median PI
brm_m5.1$criteria$loo_R2 %>% tidybayes::mode_hdi(.width = .5) # Mode HPDI

# Compare

brm_m5.2_priors <- c(prior(normal(0, 0.2), class = Intercept),
                     prior(normal(0, 0.5), class = b, coef = M),
                     prior(exponential(1), class = sigma))

brm_m5.2 <- brm(D ~ M + 1,
                data = d,
                family = gaussian,
                prior = brm_m5.2_priors,
                iter = 1000, warmup = 500, chains = 4, cores = 4, seed = 2021)

brm_m5.2 <- add_criterion(brm_m5.2, c("waic", "loo"))

brm_m5.3_priors <- c(prior(normal(0, 0.2), class = Intercept),
                     prior(normal(0, 0.5), class = b, coef = A),
                     prior(normal(0, 0.5), class = b, coef = M),
                     prior(exponential(1), class = sigma))

brm_m5.3 <- brm(D ~ A + M + 1,
                data = d,
                family = gaussian,
                prior = brm_m5.3_priors,
                iter = 1000, warmup = 500, chains = 4, cores = 4, seed = 2021)

brm_m5.3 <- add_criterion(brm_m5.3, c("waic", "loo"))

w <- loo::loo_compare(brm_m5.1, brm_m5.2, brm_m5.3, criterion = "waic")
w
print(w, simplify = F)

loo::loo_model_weights(brm_m5.1, brm_m5.2, brm_m5.3, weights  = "waic")

l <- loo::loo_compare(brm_m5.1, brm_m5.2, brm_m5.3, criterion = "loo")
l
print(l, simplify = F)

loo::loo_model_weights(brm_m5.1, brm_m5.2, brm_m5.3, weights  = "loo")
