# Ch17 - (Normal) Hierarchical Models with Predictors

# Load packages
library(bayesrules)
library(tidyverse)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(broom.mixed)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)


# Data --------------------------------------------------------------------


# Load data
data(cherry_blossom_sample)
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

# 17.1 First steps: Complete pooling --------------------------------------

complete_pooled_model <- stan_glm(
  net ~ age, 
  data = running, family = gaussian, 
  prior_intercept = normal(0, 2.5, autoscale = TRUE),
  prior = normal(0, 2.5, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

prior_summary(complete_pooled_model)

model_summary <- tidy(complete_pooled_model, 
                      conf.int = TRUE, conf.level = 0.80)
model_summary

ggplot(running, aes(age, net)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)


# 17.2 Hierarchical model with varying intercepts -------------------------


# 17.2.3 Tuning the prior -------------------------------------------------


# prior_intercept
running %>% filter(age == 55) %>% summarise(avg_net = mean(net), sd_net = sd(net))

# prior
# We’re pretty sure that the typical runner’s net time in the 10-mile race will, 
# on average, increase over time. We’re not very sure about the rate of this 
# increase, but think it’s likely between 0.5 and 4.5 minutes per year. 

# We specify prior_PD = TRUE to indicate that we wish to simulate parameters 
# from the prior, not posterior, models.
running_model_1_prior <- stan_glmer(
  net ~ age + (1 | runner), 
  data = running, family = gaussian,
  prior_intercept = normal(100, 10),
  prior = normal(2.5, 1), 
  prior_aux = exponential(1, autoscale = TRUE),
  prior_covariance = decov(reg = 1, conc = 1, shape = 1, scale = 1),
  chains = 4, iter = 5000*2, seed = 84735, 
  prior_PD = TRUE)

# we show just 4 prior plausible scenarios of what the mean regression models,  
# might look like for our 36 runners
set.seed(84735)
running %>% 
  add_fitted_draws(running_model_1_prior, n = 4) %>%
  ggplot(aes(x = age, y = net)) +
  geom_line(aes(y = .value, group = paste(runner, .draw))) + 
  facet_wrap(~ .draw)

# we also simulate 100 datasets of race outcomes from the prior model, across a 
# variety of runners and ages. 
running %>%
  add_predicted_draws(running_model_1_prior, n = 100) %>%
  ggplot(aes(x = net)) +
  geom_density(aes(x = .prediction, group = .draw)) +
  xlim(-100,300)


# 17.2.4 Posterior simulation & analysis ----------------------------------

ggplot(running, aes(age, net)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(running$runner)

# Simulate the posterior
running_model_1 <- update(running_model_1_prior, prior_PD = FALSE)

# Check the prior specifications
prior_summary(running_model_1)

# Markov chain diagnostics
mcmc_trace(running_model_1)
mcmc_dens_overlay(running_model_1)
mcmc_acf(running_model_1)
neff_ratio(running_model_1)
rhat(running_model_1)

tidy_summary_1 <- tidy(running_model_1, effects = "fixed",
                       conf.int = TRUE, conf.level = 0.80)
tidy_summary_1

# Posterior summaries of runner-specific intercepts
runner_summaries_1 <- running_model_1 %>%
  spread_draws(`(Intercept)`, b[,runner]) %>% 
  mutate(runner_intercept = `(Intercept)` + b) %>% 
  select(-`(Intercept)`, -b) %>% 
  median_qi(.width = 0.80) %>% 
  select(runner, runner_intercept, .lower, .upper)


# 17.3 Hierarchical model with varying intercepts & slopes ----------------

# Plot runner-specific models in the data
running %>% 
  filter(runner %in% c("4", "5", "20", "29")) %>% 
  ggplot(., aes(x = age, y = net)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_grid(~ runner)

ggplot(running, aes(x = age, y = net, group = runner)) + 
  geom_smooth(method = "lm", se = FALSE, size = 0.5)

# 17.3.3 Posterior simulation & analysis ----------------------------------

running_model_2 <- stan_glmer(
  net ~ age + (age | runner),
  data = running, family = gaussian,
  prior_intercept = normal(100, 10),
  prior = normal(2.5, 1), 
  prior_aux = exponential(1, autoscale = TRUE),
  prior_covariance = decov(reg = 1, conc = 1, shape = 1, scale = 1),
  chains = 4, iter = 5000*2, seed = 84735, adapt_delta = 0.99999 # <= SLOOOWW
)

# Confirm the prior model specifications
prior_summary(running_model_2)

# Quick summary of global regression parameters
tidy(running_model_2, effects = "fixed", conf.int = TRUE, conf.level = 0.80)

# spread_draws() uses b[term, runner] to grab the chains for all runner-specific 
# parameters. As usual now, these chains correspond to  b_oj and b_1j, the 
# differences between the runner-specific vs global intercepts and age coefficients.
# Get MCMC chains for the runner-specific intercepts & slopes
runner_chains_2 <- running_model_2 %>%
  spread_draws(`(Intercept)`, b[term, runner], `age`) %>% 
  pivot_wider(names_from = term, names_glue = "b_{term}",
              values_from = b) %>% 
  mutate(runner_intercept = `(Intercept)` + `b_(Intercept)`,
         runner_age = age + b_age)

runner_summaries_2 <- runner_chains_2 %>% 
  group_by(runner) %>% 
  summarize(runner_intercept = median(runner_intercept),
            runner_age = median(runner_age))

ggplot(running, aes(y = net, x = age, group = runner)) + 
  geom_abline(data = runner_summaries_2, color = "gray",
              aes(intercept = runner_intercept, slope = runner_age)) + 
  lims(x = c(50, 61), y = c(50, 135))

#  We were slightly surprised. The slopes do differ, but not as drastically as 
# we expected. But then we remembered – shrinkage! Consider sample runners 1 and 
# 10. Their posteriors suggest that, on average, runner 10’s running time 
# increases by just 1.06 minute per year, whereas runner 1’s increases by 1.75 
# minutes per year:
runner_summaries_2 %>% 
  filter(runner %in% c("runner:1", "runner:10"))


# 17.3.3.2 Posterior analysis of within- and between-group variabi --------

tidy(running_model_2, effects = "ran_pars")
