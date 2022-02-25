library(rstan)
library(bayesplot)
library(tidybayes)
library(tidyverse)

sales_ts <- readRDS(here::here("STAN", "data", "sales_ts.Rds")) %>% 
  gather(ctry_segment, value, -date) %>% 
  separate(ctry_segment, into = c("country", "segment")) %>% 
  mutate_at(vars(country, segment), ~ factor(.))

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

y <- sales_ts %>% 
  filter(country == "IT", segment == "existing") %>% 
  pull(value) %>% 
  scale() %>% 
  as.numeric()

fit_matrix_llt_seas <- stan(here::here("STAN", "stan_ts", "matrix_llt_seas.stan"), 
                            data = list(y = y,
                                        N = length(y),
                                        S = 52),
                            iter = 200, 
                            chains = 4)

plot(y, type = "l")

rstan::extract(fit_matrix_llt_seas, pars = c("y_hat_llt", "x_seas")) %>% 
  map_df(colMeans) %>% 
  mutate(y_hat = y_hat_llt + x_seas) %>% 
  pull(y_hat) %>% 
  lines(col = "blue") 

as.data.frame(fit_matrix_llt_seas) %>% 
  select(starts_with("mu_y[")) %>% 
  colMeans() %>% 
  lines(col = "blue")

## Scaled intetrnally
y <- sales_ts %>% 
  filter(country == "IT", segment == "existing") %>% 
  pull(value) 

fit_matrix_llt_seas <- stan(here::here("STAN", "stan_ts", "matrix_scaled_llt_seas.stan"), 
                            data = list(y = y,
                                        N = length(y),
                                        S = 52),
                            iter = 2000, 
                            chains = 4)

plot(y, type = "l")

rstan::extract(fit_matrix_llt_seas, pars = c("mu")) %>% 
  map_df(colMeans) %>% 
  pull(mu) %>% 
  lines(col = "blue") 

as.data.frame(fit_matrix_llt_seas) %>% 
  select(starts_with("mu_y[")) %>% 
  colMeans() %>% 
  lines(col = "blue")
