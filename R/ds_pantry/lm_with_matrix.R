library(tidyverse)
library(bayesplot)
library(tidybayes)
library(broom.mixed)
library(rstan)
library(rethinking)

data(Howell1)

d <- Howell1 %>% as_tibble()
d2 <- Howell1 %>% filter(age >= 18) %>% as_tibble()


options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


# Modelo básico -----------------------------------------------------------

stan_m1 <- stan_model(file = "R/ds_pantry/lm_with_matrix.stan")

fit_stan_m1 <- sampling(stan_m1, 
                        iter = 10000, chains = 1, 
                        # control = list(adapt_delta = 0.99),
                        data = list(y = d$height,
                                    X = model.matrix(height ~ ., d),
                                    n = nrow(d),
                                    m = 4),
                        verbose = TRUE)

print(fit_stan_m1, pars = c("betas"))

fit_lm <- lm(height ~ ., d)
summary(fit_lm)
plot(d$height, fitted(fit_lm))

plot(d$height, extract(fit_stan_m1, "mu")[[1]] %>% apply(2, mean))

# Escalado de datos -------------------------------------------------------

stan_m2 <- stan_model(file = "R/ds_pantry/lm_with_matrix2.stan")

fit_stan_m2 <- sampling(stan_m2, 
                        iter = 10000, chains = 1, 
                        # control = list(adapt_delta = 0.99),
                        data = c(compose_data(d), intcpt = TRUE),
                        verbose = TRUE)

print(fit_stan_m2, pars = c("beta_intcpt", "beta_weight", "beta_age", "beta_male"))
traceplot(fit_stan_m2, pars = c("beta_intcpt", "beta_weight", "beta_age", "beta_male"))

fit_lm <- lm(height ~ ., d)
summary(fit_lm)

plot(d$height, fitted(fit_lm))

plot(d$height, extract(fit_stan_m2, "mu")[[1]] %>% apply(2, mean))

plot(fitted(fit_lm), extract(fit_stan_m2, "mu")[[1]] %>% apply(2, mean))


# Tidybayes

tidybayes::spread_draws(fit_stan_m2, beta_intcpt, beta_weight, beta_age, beta_male)

tidybayes::gather_draws(fit_stan_m2, beta_intcpt, beta_weight, beta_age, beta_male)

tidybayes::gather_draws(fit_stan_m2, beta_intcpt, beta_weight, beta_age, beta_male) %>% 
  ggplot(aes(x=.value)) + 
  geom_density() + 
  facet_wrap(~.variable, scales = "free")

tidybayes::median_qi(tidybayes::spread_draws(fit_stan_m2, beta_intcpt, beta_weight, beta_age, beta_male), .width = .89)

add_epred_draws(newdata = d, object = fit_stan_m2)
epred_draws(object = fit_stan_m2, newdata = d)

add_linpred_draws(newdata = d, object = fit_stan_m2)
linpred_draws(object = fit_stan_m2, newdata = d)

add_predicted_draws(newdata = d, object = fit_stan_m2)
predicted_draws(object = fit_stan_m2, newdata = d)

add_residual_draws(newdata = d, object = fit_stan_m2)
residual_draws(object = fit_stan_m2, newdata = d)

rstantools::posterior_epred(fit_stan_m2)


# Reparametrización QR ----------------------------------------------------

x <- model.matrix(height ~.-1 , d)

stan_m3 <- stan_model(file = "R/ds_pantry/lm_with_matrix3.stan")

fit_stan_m3 <- sampling(stan_m3, 
                        iter = 10000, chains = 1, 
                        # control = list(adapt_delta = 0.99),
                        data = c(compose_data(x),
                                 y = list(d$height),
                                 N = nrow(x), K = ncol(x)),
                        verbose = TRUE)

print(fit_stan_m3, pars = c("beta"))
traceplot(fit_stan_m3, pars = c("beta_intcpt", "beta_weight", "beta_age", "beta_male"))

fit_lm <- lm(height ~ ., d)
summary(fit_lm)

plot(d$height, fitted(fit_lm))

plot(d$height, extract(fit_stan_m2, "mu")[[1]] %>% apply(2, mean))

plot(fitted(fit_lm), extract(fit_stan_m3, "mu")[[1]] %>% apply(2, mean))


###

d %>% group_by(male) %>% 
  summarise_all(list(avg = ~ mean(.x), sd = ~ sd(.x)),
            .names = "{.col}_{.fn}")
