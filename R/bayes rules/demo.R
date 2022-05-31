library(tidyverse)

N_G <- 4

# set.seed(9999)
set.seed(9999)
a_0 <- rnorm(1, 5, 10)
a_1 <- rnorm(1, 0, 2)
b_0 <- rnorm(1, 0, 3)
b_1 <- rnorm(1, 0, 5)


conf_grupos <- tibble(G = 1:N_G, 
                      beta_0 = rnorm(N_G, a_0, 10),
                      beta_1 = rnorm(N_G, b_0, 2),
                      sigma_y = 50)

X <- seq(0, 10, by = .1)

n_j <- sample.int(length(X)-3, N_G) + 3

complete_dataset <- lapply(seq_along(n_j), function(i) {
  X = sort(sample(X, n_j[i]))
  tibble(G = i, I = 1:n_j[i], X)
}) %>% 
  bind_rows() %>% 
  left_join(conf_grupos, by = "G") %>% 
  mutate(y = rnorm(length(X), beta_0 + beta_1 * X, sigma_y),
         G = factor(G))

dataset <- complete_dataset %>% 
  select(G, I, X, y)

ggplot(data = dataset) +
  geom_density(mapping = aes(x = y, fill = G, alpha = .1))

p <- ggplot() +
  geom_point(data = dataset, mapping = aes(x = X, y = y, group = G, color = G)) 
p

p <- p + geom_smooth(data = dataset, aes(x = X, y = y), 
                     color = "black", size = 2, method = "lm", se = FALSE) 
p

p <- p + geom_smooth(data = dataset, mapping = aes(x = X, y = y, group = G, color = G), 
                     method = "lm", se = FALSE)

p
  
library(lme4)

me_lm <- lmer(y ~ X + (X | G), dataset)
summary(me_lm)

coef_me_lm <- coef(me_lm)$G %>% 
  rownames_to_column("G") %>% 
  rename(intercept = `(Intercept)`, slope = X) #%>% 
    # gather(coef, value, -G)

p +  
  # scale_x_continuous(name="x", limits=c(0,75)) +
  # scale_y_continuous(name="y", limits=c(-50,500)) +
geom_abline(data = coef_me_lm, 
            mapping = aes(intercept = intercept, slope = slope, color = G), 
            linetype = 2, size = 1.2)


library(rstanarm)
library(tidybayes)
library(broom.mixed)
library(bayes_rules)
library(bayesplot)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

bme_lm <- stan_glmer(y ~ X + (X | G), 
                     family =  gaussian, data = dataset, 
                     chains = 4, iter = 5000*2, seed = 84735)

# Prior specification
prior_summary(bme_lm)

# MCMC diagnostics
mcmc_trace(bme_lm, size = 0.1)
mcmc_dens_overlay(bme_lm)
mcmc_acf(bme_lm)
neff_ratio(bme_lm)
rhat(bme_lm)

pp_check(bme_lm) 

bme_lm_df <- as.data.frame(bme_lm)

G_chains_intercept <- bme_lm %>%
  spread_draws(`(Intercept)`, b[,G]) %>% 
  mutate(mu_j = `(Intercept)` + b) 

G_chains_intercept_scaled <- G_chains_intercept %>% 
  select(-`(Intercept)`, -b) %>% 
  mean_qi(.width = 0.80) %>% 
  mutate(artist = fct_reorder(G, mu_j))

G_chains_X <- bme_lm %>%
  spread_draws(X, b[,G]) %>% 
  mutate(mu_j = X + b) 

G_chains_X_scaled <- G_chains_X %>% 
  select(-X, -b) %>% 
  mean_qi(.width = 0.80) %>% 
  mutate(artist = fct_reorder(G, mu_j))

broom.mixed::tidy(bme_lm, effects = "fixed", conf.int = TRUE, conf.level = .8)

broom.mixed::tidy(bme_lm, effects = "ran_pars", conf.int = TRUE, conf.level = .8)

bayesrules::prediction_summary(model = bme_lm, data = dataset)
loo(bme_lm)

####
# Partially pooled
me_lm <- lmer(y ~ X + (X | G), dataset)
summary(me_lm)
coef(me_lm)

# Pooled
me_lm_pooled <- lm(y ~ X , dataset)
summary(me_lm_pooled)
coef(me_lm_pooled)

# No pooling
me_lm_no_pooling <- lm(y ~ X + G + X*G, dataset)
# lm(y ~ X:G + G, dataset)
summary(me_lm_no_pooling)
coef(me_lm_no_pooling)
