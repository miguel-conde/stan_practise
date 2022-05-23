# Ch 16 - (Normal) Hierarchical Models without Predictors

library(tidyverse)
library(bayesrules)

library(rstanarm)
library(bayesplot)
library(tidybayes)
library(broom.mixed)
library(forcats)

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

# Data --------------------------------------------------------------------


data(spotify)

spotify

glimpse(spotify)

spotify <- spotify %>% 
  select(artist, title, popularity) %>% 
  mutate(artist = fct_reorder(artist, popularity, .fun = "mean"))

spotify
glimpse(spotify)
nlevels(spotify$artist)

artist_means <- spotify %>% 
  group_by(artist) %>% 
  summarise(count = n(), popularity = mean(popularity))

artist_means %>% slice(1:2, 43:44)


# Complete pooled model ---------------------------------------------------

spotify %>% 
  ggplot(aes(x = popularity)) +
  geom_density()

# j -> Artist
# i -> i-th song of artist i
#
# Y_i_j | mu, sigma ~ N(mu, sigma^2)
#                mu ~ N(50, )
#             sigma ~ Exp()

# Solo tenemos un a priori (mu aprox 50), así que los otros se dejan con sus
# valores por defecto (weakly informative)

spotify_complete_pooled <- 
  stan_glm(formula = popularity ~ 1,
           data = spotify,
           family = gaussian,
           prior_intercept = normal(50, 2.5, autoscale = TRUE),
           prior_aux = exponential(1, autoscale = TRUE),
           chains = 4, iter = 5000*2, seed = 84735)

prior_summary(spotify_complete_pooled)

complete_summary <- tidy(spotify_complete_pooled,
                         effects = c("fixed", "aux"),
                         conf.int = TRUE, conf.level = .80)
complete_summary

set.seed(84735)
prediction_complete <- posterior_predict(spotify_complete_pooled, newdata = artist_means)

ppc_intervals(y = artist_means$popularity, yrep = prediction_complete, prob_outer = .80) +
  ggplot2::scale_x_continuous(labels = artist_means$artist, 
                              breaks = 1:nrow(artist_means)) +
  xaxis_text(angle = 90, hjust = 1)


# No pooled model ---------------------------------------------------------

ggplot(spotify, aes(x = popularity, group = artist)) +
  geom_density()

# Y_i_j | mu_j, sigma ~ N(mu_j, sigma^2)
#                mu_j ~ N(50, )
#             sigma ~ Exp()

spotify_no_pooled <- 
  stan_glm(formula = popularity ~ artist - 1,
           data = spotify,
           family = gaussian,
           prior_intercept = normal(50, 2.5, autoscale = TRUE),
           prior_aux = exponential(1, autoscale = TRUE),
           chains = 4, iter = 5000*2, seed = 84735)

prior_summary(spotify_no_pooled)

complete_summary <- tidy(spotify_no_pooled,
                         effects = c("fixed", "aux"),
                         conf.int = TRUE, conf.level = .80)
complete_summary

set.seed(84735)
prediction_complete <- posterior_predict(spotify_no_pooled, newdata = artist_means)

ppc_intervals(y = artist_means$popularity, yrep = prediction_complete, prob_outer = .80) +
  ggplot2::scale_x_continuous(labels = artist_means$artist, 
                              breaks = 1:nrow(artist_means)) +
  xaxis_text(angle = 90, hjust = 1)

# Hierarchical ------------------------------------------------------------

## Posterior Simulation

spotify_hierarchical <- 
  stan_glmer(formula = popularity ~ 1 + (1 | artist),
             data = spotify,
             family = gaussian,
             prior_intercept = normal(50, 2.5, autoscale = TRUE),
             prior_aux = exponential(1, autoscale = TRUE),
             prior_covariance = decov(reg = 1, conc = 1, shape = 1, scale = 1),
             chains = 4, iter = 5000*2, seed = 84735)

spotify_hierarchical

prior_summary(spotify_hierarchical)

mcmc_trace(spotify_hierarchical)
mcmc_dens_overlay(spotify_hierarchical)
mcmc_acf(spotify_hierarchical)
neff_ratio(spotify_hierarchical)
rhat(spotify_hierarchical)

pp_check(spotify_hierarchical) + xlab("popularity")

spotify_hierarchical_df <- spotify_hierarchical %>% as.data.frame()
spotify_hierarchical_df %>% names()

## Posterior analysis of global parameters

complete_summary <- tidy(spotify_hierarchical,
                         effects = c("fixed", "aux", "ran_pars"),
                         conf.int = TRUE, conf.level = .80)
complete_summary

tidy(spotify_hierarchical,
     effects = c("fixed"),
     conf.int = TRUE, conf.level = .80)

tidy(spotify_hierarchical,
     effects = c("ran_pars"),
     conf.int = TRUE, conf.level = .80)

# The popularity levels among multiple songs by the same
# artist tend to have a moderate correlation near 0.54:
15.2^2 / (15.2^2 + 14^2)

# Thinking of this another way, 54% of the variability in song popularity is
# explained by differences between artists, whereas 46% is explained by
# differences among the songs within each artist:
14 ^2 / (15.2^2 + 14^2)

set.seed(84735)
prediction_complete <- posterior_predict(spotify_hierarchical, newdata = artist_means)

ppc_intervals(y = artist_means$popularity, yrep = prediction_complete, prob_outer = .80) +
  ggplot2::scale_x_continuous(labels = artist_means$artist, 
                              breaks = 1:nrow(artist_means)) +
  xaxis_text(angle = 90, hjust = 1)

## Posterior analysis of group-specifica parameters

artists_summary <- tidy(spotify_hierarchical,
                        effects = c("ran_vals"),
                        conf.int = TRUE, conf.level = .80)
artists_summary # Son relativos a la valoracion del artista medio

artist_chains <- spotify_hierarchical %>% 
  spread_draws(`(Intercept)`, b[, artist]) %>% 
  mutate(mu_j = `(Intercept)` + b)

artists_summary_scaled <- artist_chains %>% 
  select(artist, mu_j) %>% 
  mean_qi(.width = 0.8) %>% 
  mutate(artist = fct_reorder(artist, mu_j))

ggplot(artists_summary_scaled,
       aes(x = artist, y = mu_j, ymin = .lower, ymax = .upper)) +
  geom_pointrange() +
  xaxis_text(angle = 90, hjust = 1)

## Posterior prediction

# First consider the posterior prediction for an observed group or artist

set.seed(84735)
ocean_pps <- spotify_hierarchical_df %>% 
  select(`(Intercept)`, `b[(Intercept) artist:Frank_Ocean]`, sigma) %>% 
  mutate(mu_ocean = `(Intercept)` + `b[(Intercept) artist:Frank_Ocean]`,
         .keep = "unused") %>% 
  mutate(y_ocean = rnorm(20000, mu_ocean, sigma))

ocean_pps %>% mean_qi(y_ocean, .width = .8)

# Naturally, we can be much more certain about Ocean's underlying mean song 
# popularity than in the popularity of any single Ocean song.

artists_summary_scaled %>% 
  filter(str_detect(artist, "Ocean"))

# Next consider posterior prediction for a yet unobserved group

